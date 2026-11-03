from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from math import cos, pi, sin
from random import Random
from statistics import mean
from typing import Dict, List, Tuple


# ----------------------------
# Data models
# ----------------------------


@dataclass
class Bus:
    idx: int
    name: str
    base_kv: float
    load_mw_profile: List[float]
    renewable_sites: List[str] = field(default_factory=list)


@dataclass
class Line:
    name: str
    from_bus: int
    to_bus: int
    x_pu: float
    thermal_limit_mw: float
    in_service: bool = True


@dataclass
class RenewablePlant:
    name: str
    plant_type: str
    bus: int
    p_max_mw: float
    profile: List[float]
    inertia_emulation_mw_s: float = 0.0


@dataclass
class BatteryStorage:
    name: str
    bus: int
    p_charge_max_mw: float
    p_discharge_max_mw: float
    e_max_mwh: float
    soc_mwh: float
    eta_charge: float = 0.95
    eta_discharge: float = 0.95


@dataclass
class SimulationConfig:
    base_mva: float = 100.0
    hours: int = 24
    dt_s: float = 0.1
    transient_window_s: float = 25.0
    nominal_frequency_hz: float = 50.0
    # Stress/fault controls for high-renewable edge-case studies.
    fault_hours: Tuple[int, ...] = (8, 13, 19, 21)
    grid_forming_reserve_mw: float = 1750.0
    emergency_shed_fraction: float = 0.12
    inverter_trip_loading_pu: float = 1.35
    cascade_trip_loading_pu: float = 1.55
    blackout_loading_threshold_pu: float = 2.30
    hydrogen_flexible_load_max_mw: float = 150.0
    hydrogen_bus: int = 5


# ----------------------------
# Utility math
# ----------------------------


def gaussian_elimination(a: List[List[float]], b: List[float]) -> List[float]:
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(m[r][col]))
        if abs(m[pivot][col]) < 1e-12:
            raise ValueError("Singular matrix encountered in load flow solve.")
        m[col], m[pivot] = m[pivot], m[col]

        pivot_val = m[col][col]
        for j in range(col, n + 1):
            m[col][j] /= pivot_val

        for r in range(n):
            if r == col:
                continue
            factor = m[r][col]
            for j in range(col, n + 1):
                m[r][j] -= factor * m[col][j]

    return [m[i][n] for i in range(n)]


def clamp(v: float, vmin: float, vmax: float) -> float:
    return max(vmin, min(vmax, v))


# ----------------------------
# Core system model
# ----------------------------


class RenewableGridModel:
    def __init__(self, config: SimulationConfig, seed: int = 17):
        self.cfg = config
        self.rng = Random(seed)

        self.buses: List[Bus] = self._build_buses()
        self.lines: List[Line] = self._build_lines()
        self.plants: List[RenewablePlant] = self._build_plants()
        self.storage: List[BatteryStorage] = self._build_storage()
        self.fault_events: Dict[int, Dict[str, object]] = self._build_fault_events()

        self.bus_count = len(self.buses)
        self.hourly_results: List[Dict[str, float]] = []

    def _daily_load_shape(self) -> List[float]:
        # A realistic-ish normalized daily shape with evening peak.
        shape = []
        for h in range(self.cfg.hours):
            morning = 0.15 * sin((h - 6) * pi / 12)
            evening = 0.25 * sin((h - 16) * pi / 8)
            base = 0.85 + max(0.0, morning) + max(0.0, evening)
            shape.append(base)
        peak = max(shape)
        return [s / peak for s in shape]

    def _wind_profile(self) -> List[float]:
        out = []
        for h in range(self.cfg.hours):
            trend = 0.45 + 0.2 * sin((h + 3) * pi / 12)
            noise = self.rng.uniform(-0.08, 0.08)
            out.append(clamp(trend + noise, 0.05, 0.95))
        return out

    def _solar_profile(self) -> List[float]:
        out = []
        for h in range(self.cfg.hours):
            x = (h - 6) / 12
            day = max(0.0, sin(pi * x))
            cloud = 1.0 - self.rng.uniform(0.0, 0.25)
            out.append(clamp(day * cloud, 0.0, 1.0))
        return out

    def _firm_renewable_profile(self, cf: float) -> List[float]:
        return [cf for _ in range(self.cfg.hours)]

    def _build_buses(self) -> List[Bus]:
        shape = self._daily_load_shape()
        bases = [
            (0, "NorthHub", 220.0, 380.0),
            (1, "CoastalWind", 220.0, 260.0),
            (2, "SolarBasin", 132.0, 210.0),
            (3, "HydroValley", 220.0, 310.0),
            (4, "UrbanCenter", 400.0, 620.0),
            (5, "IndustrialEast", 220.0, 470.0),
            (6, "GeoSouth", 132.0, 190.0),
            (7, "Intertie", 400.0, 280.0),
        ]

        buses = []
        for idx, name, kv, peak_mw in bases:
            load = [peak_mw * p for p in shape]
            buses.append(Bus(idx=idx, name=name, base_kv=kv, load_mw_profile=load))
        return buses

    def _build_lines(self) -> List[Line]:
        return [
            Line("L01", 0, 1, x_pu=0.18, thermal_limit_mw=260),
            Line("L12", 1, 2, x_pu=0.15, thermal_limit_mw=240),
            Line("L23", 2, 3, x_pu=0.14, thermal_limit_mw=240),
            Line("L34", 3, 4, x_pu=0.10, thermal_limit_mw=360),
            Line("L45", 4, 5, x_pu=0.11, thermal_limit_mw=350),
            Line("L56", 5, 6, x_pu=0.16, thermal_limit_mw=220),
            Line("L67", 6, 7, x_pu=0.12, thermal_limit_mw=250),
            Line("L70", 7, 0, x_pu=0.09, thermal_limit_mw=320),
            Line("L04", 0, 4, x_pu=0.20, thermal_limit_mw=200),
            Line("L17", 1, 7, x_pu=0.13, thermal_limit_mw=230),
            Line("L25", 2, 5, x_pu=0.17, thermal_limit_mw=220),
            Line("L36", 3, 6, x_pu=0.14, thermal_limit_mw=240),
        ]

    def _build_plants(self) -> List[RenewablePlant]:
        wind = self._wind_profile()
        solar = self._solar_profile()
        hydro = self._firm_renewable_profile(0.72)
        bio = self._firm_renewable_profile(0.78)
        geo = self._firm_renewable_profile(0.86)

        plants = [
            RenewablePlant("Wind_NorthSea", "wind", 1, 280.0, wind, inertia_emulation_mw_s=180),
            RenewablePlant("Wind_CliffEdge", "wind", 7, 160.0, wind, inertia_emulation_mw_s=105),
            RenewablePlant("Solar_Delta", "solar", 2, 300.0, solar, inertia_emulation_mw_s=0),
            RenewablePlant("Solar_Desert", "solar", 5, 250.0, solar, inertia_emulation_mw_s=0),
            RenewablePlant("Hydro_Main", "hydro", 3, 340.0, hydro, inertia_emulation_mw_s=420),
            RenewablePlant("Biomass_CHP", "biomass", 0, 120.0, bio, inertia_emulation_mw_s=90),
            RenewablePlant("Geothermal_Core", "geothermal", 6, 190.0, geo, inertia_emulation_mw_s=210),
        ]

        for p in plants:
            self.buses[p.bus].renewable_sites.append(p.name)
        return plants

    def _build_storage(self) -> List[BatteryStorage]:
        return [
            BatteryStorage(
                name="BESS_Urban",
                bus=4,
                p_charge_max_mw=150,
                p_discharge_max_mw=180,
                e_max_mwh=600,
                soc_mwh=290,
            ),
            BatteryStorage(
                name="BESS_WindPort",
                bus=1,
                p_charge_max_mw=90,
                p_discharge_max_mw=110,
                e_max_mwh=320,
                soc_mwh=170,
            ),
        ]

    def _build_fault_events(self) -> Dict[int, Dict[str, object]]:
        events: Dict[int, Dict[str, object]] = {}
        fault_lines = ["L34", "L25", "L70", "L45", "L12", "L67"]
        fault_types = ["3PH", "LL", "SLG", "3PH"]

        for idx, hour in enumerate(self.cfg.fault_hours):
            line_name = fault_lines[idx % len(fault_lines)]
            ftype = fault_types[idx % len(fault_types)]
            events[int(hour)] = {
                "type": ftype,
                "line_outages": [line_name],
                "load_spike": 0.04 + 0.02 * (idx % 3),
                "wind_derate": 0.16 + 0.04 * (idx % 2),
                "solar_derate": 0.22 + 0.05 * (idx % 2),
                "voltage_sag": 0.08 + 0.03 * (idx % 2),
            }
        return events

    def _line_service_map(self) -> Dict[str, bool]:
        return {ln.name: ln.in_service for ln in self.lines}

    def _lines_from_service(self, service_map: Dict[str, bool]) -> List[Line]:
        return [
            Line(
                name=ln.name,
                from_bus=ln.from_bus,
                to_bus=ln.to_bus,
                x_pu=ln.x_pu,
                thermal_limit_mw=ln.thermal_limit_mw,
                in_service=service_map.get(ln.name, False),
            )
            for ln in self.lines
        ]

    def _incident_line_loading(
        self,
        bus: int,
        line_flows: Dict[str, float],
        service_map: Dict[str, bool],
    ) -> float:
        loading = 0.0
        for ln in self.lines:
            if not service_map.get(ln.name, False):
                continue
            if ln.from_bus == bus or ln.to_bus == bus:
                loading = max(loading, abs(line_flows.get(ln.name, 0.0)) / ln.thermal_limit_mw)
        return loading

    def _simulate_hour_edge_case(
        self,
        hour: int,
        soc_by_name: Dict[str, float],
        line_service: Dict[str, bool],
    ) -> Dict[str, object]:
        event = self.fault_events.get(hour)
        plant_mw = self._plant_output(hour)

        if event:
            for p in self.plants:
                if p.plant_type == "wind":
                    plant_mw[p.name] *= (1.0 - float(event["wind_derate"]))
                if p.plant_type == "solar":
                    plant_mw[p.name] *= (1.0 - float(event["solar_derate"]))
            for ln_name in event["line_outages"]:
                line_service[ln_name] = False

        total_load = sum(b.load_mw_profile[hour] for b in self.buses)
        if event:
            total_load *= (1.0 + float(event["load_spike"]))

        # Hydrogen electrolyzer as flexible demand: absorbs excess renewable, backs off under stress.
        total_renew = sum(plant_mw.values())
        mismatch_pre_h2 = total_load - total_renew
        h2_load = 0.0
        if mismatch_pre_h2 < -80.0:
            h2_load = min(self.cfg.hydrogen_flexible_load_max_mw, -0.55 * mismatch_pre_h2)
        elif mismatch_pre_h2 > 120.0:
            h2_load = 0.0
        else:
            h2_load = 0.22 * self.cfg.hydrogen_flexible_load_max_mw
        total_load += h2_load

        mismatch = total_load - total_renew
        storage_dispatch: Dict[str, float] = {s.name: 0.0 for s in self.storage}
        for s in self.storage:
            soc = soc_by_name[s.name]
            if mismatch > 1e-6:
                max_energy_limited = soc * s.eta_discharge
                p = min(s.p_discharge_max_mw, mismatch, max_energy_limited)
                soc_by_name[s.name] = soc - p / s.eta_discharge
                storage_dispatch[s.name] = p
                mismatch -= p
            elif mismatch < -1e-6:
                room = s.e_max_mwh - soc
                p = min(s.p_charge_max_mw, -mismatch, room / s.eta_charge)
                soc_by_name[s.name] = soc + p * s.eta_charge
                storage_dispatch[s.name] = -p
                mismatch += p

        storage_net = sum(storage_dispatch.values())
        mismatch_after_storage = total_load - (total_renew + storage_net)
        dr = self._demand_response(mismatch_after_storage, total_load)
        mismatch_after_dr = mismatch_after_storage - dr

        emergency_shed = 0.0
        if event and mismatch_after_dr > 0.0:
            emergency_shed = min(
                self.cfg.emergency_shed_fraction * total_load,
                0.60 * mismatch_after_dr,
            )
            mismatch_after_dr -= emergency_shed

        curtailment = max(0.0, -mismatch_after_dr)
        mismatch_after_dr = max(0.0, mismatch_after_dr)
        balancer = min(self.cfg.grid_forming_reserve_mw, mismatch_after_dr)
        unserved = max(0.0, mismatch_after_dr - balancer)

        bus_generation = [0.0 for _ in range(self.bus_count)]
        for p in self.plants:
            bus_generation[p.bus] += plant_mw[p.name]
        for s in self.storage:
            bus_generation[s.bus] += storage_dispatch[s.name]
        bus_generation[0] += balancer

        p_inj = []
        for i, b in enumerate(self.buses):
            effective_load = b.load_mw_profile[hour] * (1.0 + float(event["load_spike"])) if event else b.load_mw_profile[hour]
            if i == self.cfg.hydrogen_bus:
                effective_load += h2_load
            p_inj.append(bus_generation[i] - effective_load + (dr + emergency_shed if i == 4 else 0.0))

        cascaded_lines: List[str] = []
        islanded = False
        for _ in range(3):
            case_lines = self._lines_from_service(line_service)
            try:
                _, line_flows = self.dc_load_flow(p_inj, case_lines)
            except ValueError:
                islanded = True
                line_flows = {ln.name: 0.0 for ln in self.lines}
                break

            worst_line = ""
            worst_loading = 0.0
            for ln in self.lines:
                if not line_service.get(ln.name, False):
                    continue
                loading = abs(line_flows[ln.name]) / ln.thermal_limit_mw
                if loading > worst_loading:
                    worst_loading = loading
                    worst_line = ln.name

            trigger_threshold = self.cfg.cascade_trip_loading_pu if event else 9.9
            if worst_loading <= trigger_threshold:
                break
            if worst_line:
                line_service[worst_line] = False
                cascaded_lines.append(worst_line)

        # Inverter ride-through violation and protection trips.
        inverter_tripped_mw = 0.0
        if not islanded:
            case_lines = self._lines_from_service(line_service)
            try:
                _, line_flows = self.dc_load_flow(p_inj, case_lines)
            except ValueError:
                islanded = True
                line_flows = {ln.name: 0.0 for ln in self.lines}
            for p in self.plants:
                if p.plant_type not in {"wind", "solar"}:
                    continue
                local_loading = self._incident_line_loading(p.bus, line_flows, line_service)
                sag = float(event["voltage_sag"]) if event else 0.0
                if local_loading > self.cfg.inverter_trip_loading_pu or sag > 0.09:
                    cut = plant_mw[p.name] * min(0.72, 0.25 + 0.35 * max(0.0, local_loading - 1.0) + 0.8 * sag)
                    plant_mw[p.name] -= cut
                    inverter_tripped_mw += cut

            if inverter_tripped_mw > 0.0:
                bus_generation = [0.0 for _ in range(self.bus_count)]
                for p in self.plants:
                    bus_generation[p.bus] += plant_mw[p.name]
                for s in self.storage:
                    bus_generation[s.bus] += storage_dispatch[s.name]
                bus_generation[0] += balancer
                p_inj = []
                for i, b in enumerate(self.buses):
                    effective_load = b.load_mw_profile[hour] * (1.0 + float(event["load_spike"])) if event else b.load_mw_profile[hour]
                    if i == self.cfg.hydrogen_bus:
                        effective_load += h2_load
                    p_inj.append(bus_generation[i] - effective_load + (dr + emergency_shed if i == 4 else 0.0))
                case_lines = self._lines_from_service(line_service)
                try:
                    _, line_flows = self.dc_load_flow(p_inj, case_lines)
                except ValueError:
                    islanded = True
                    line_flows = {ln.name: 0.0 for ln in self.lines}
        else:
            line_flows = {ln.name: 0.0 for ln in self.lines}

        total_renew_post_trip = sum(plant_mw.values())
        if inverter_tripped_mw > 0:
            unserved += max(0.0, inverter_tripped_mw - 0.35 * self.cfg.grid_forming_reserve_mw)

        max_loading = max((abs(line_flows[ln.name]) / ln.thermal_limit_mw for ln in self.lines), default=0.0)
        overloaded = {
            ln.name: abs(line_flows[ln.name]) / ln.thermal_limit_mw
            for ln in self.lines
            if line_service.get(ln.name, False) and abs(line_flows[ln.name]) > ln.thermal_limit_mw
        }

        voltage_index = self._voltage_proxy(p_inj) - (float(event["voltage_sag"]) if event else 0.0)
        voltage_index = clamp(voltage_index - 0.06 * (len(cascaded_lines)), 0.50, 1.05)

        inertia = self._system_inertia(hour)
        rocof_proxy = (unserved / max(inertia, 100.0)) * 12.0
        freq_proxy = self.cfg.nominal_frequency_hz - 0.012 * unserved - 0.18 * max(0.0, rocof_proxy - 0.8)
        freq_proxy = clamp(freq_proxy, 44.5, self.cfg.nominal_frequency_hz)

        blackout = bool(
            islanded
            or (
                event is not None
                and (
                    max_loading > self.cfg.blackout_loading_threshold_pu
                    or freq_proxy < 48.6
                    or voltage_index < 0.62
                    or unserved > 140.0
                )
            )
        )

        return {
            "hour": hour,
            "fault_active": float(event is not None),
            "fault_type": event["type"] if event else "NONE",
            "fault_lines": ",".join(event["line_outages"]) if event else "",
            "line_service": dict(line_service),
            "cascaded_lines": ",".join(cascaded_lines),
            "islanded": float(islanded),
            "blackout": float(blackout),
            "p_inj": p_inj,
            "line_flows": line_flows,
            "max_loading_pu": max_loading,
            "load_mw": total_load,
            "renewable_mw": total_renew_post_trip,
            "storage_net_mw": storage_net,
            "demand_response_mw": dr,
            "emergency_shed_mw": emergency_shed,
            "curtailment_mw": curtailment,
            "balancer_mw": balancer,
            "unserved_mw": unserved,
            "hydrogen_load_mw": h2_load,
            "overloaded_lines": float(len(overloaded)),
            "voltage_stability_idx": voltage_index,
            "inertia_mw_s": inertia,
            "inverter_tripped_mw": inverter_tripped_mw,
            "frequency_proxy_hz": freq_proxy,
            "rocof_proxy_hz_s": rocof_proxy,
        }

    def _plant_output(self, hour: int) -> Dict[str, float]:
        return {p.name: p.p_max_mw * p.profile[hour] for p in self.plants}

    def _dispatch_storage(self, hour: int, mismatch_mw: float) -> Dict[str, float]:
        # Positive mismatch means deficit -> discharge storage.
        dispatch: Dict[str, float] = {s.name: 0.0 for s in self.storage}

        for s in self.storage:
            if mismatch_mw > 1e-6:
                max_energy_limited = s.soc_mwh * s.eta_discharge
                p = min(s.p_discharge_max_mw, mismatch_mw, max_energy_limited)
                s.soc_mwh -= p / s.eta_discharge
                dispatch[s.name] = p
                mismatch_mw -= p
            elif mismatch_mw < -1e-6:
                room = s.e_max_mwh - s.soc_mwh
                p = min(s.p_charge_max_mw, -mismatch_mw, room / s.eta_charge)
                s.soc_mwh += p * s.eta_charge
                dispatch[s.name] = -p
                mismatch_mw += p

        return dispatch

    def _demand_response(self, deficit_mw: float, total_load_mw: float) -> float:
        # Up to 5% flexible curtailment under stress.
        max_dr = 0.05 * total_load_mw
        return min(max(deficit_mw, 0.0), max_dr)

    def _build_b_matrix(self, lines: List[Line]) -> List[List[float]]:
        n = self.bus_count
        b = [[0.0 for _ in range(n)] for _ in range(n)]

        for ln in lines:
            if not ln.in_service:
                continue
            bij = 1.0 / ln.x_pu
            i, j = ln.from_bus, ln.to_bus
            b[i][i] += bij
            b[j][j] += bij
            b[i][j] -= bij
            b[j][i] -= bij
        return b

    def dc_load_flow(self, p_inj_mw: List[float], lines: List[Line]) -> Tuple[List[float], Dict[str, float]]:
        slack = 0
        n = self.bus_count
        b = self._build_b_matrix(lines)

        reduced = []
        rhs = []
        for i in range(n):
            if i == slack:
                continue
            row = []
            for j in range(n):
                if j == slack:
                    continue
                row.append(b[i][j])
            reduced.append(row)
            rhs.append(p_inj_mw[i] / self.cfg.base_mva)

        theta_red = gaussian_elimination(reduced, rhs)
        theta = [0.0] * n
        k = 0
        for i in range(n):
            if i == slack:
                continue
            theta[i] = theta_red[k]
            k += 1

        line_flows = {}
        for ln in lines:
            if not ln.in_service:
                line_flows[ln.name] = 0.0
                continue
            f = (theta[ln.from_bus] - theta[ln.to_bus]) / ln.x_pu * self.cfg.base_mva
            line_flows[ln.name] = f

        return theta, line_flows

    def _hourly_power_balance(self, hour: int) -> Dict[str, float]:
        local_soc = {s.name: s.soc_mwh for s in self.storage}
        local_lines = self._line_service_map()
        st = self._simulate_hour_edge_case(hour, local_soc, local_lines)
        return {
            "hour": float(st["hour"]),
            "load_mw": float(st["load_mw"]),
            "renewable_mw": float(st["renewable_mw"]),
            "storage_net_mw": float(st["storage_net_mw"]),
            "demand_response_mw": float(st["demand_response_mw"]),
            "emergency_shed_mw": float(st["emergency_shed_mw"]),
            "curtailment_mw": float(st["curtailment_mw"]),
            "balancer_mw": float(st["balancer_mw"]),
            "unserved_mw": float(st["unserved_mw"]),
            "hydrogen_load_mw": float(st["hydrogen_load_mw"]),
            "fault_active": float(st["fault_active"]),
            "blackout": float(st["blackout"]),
            "overloaded_lines": float(st["overloaded_lines"]),
            "max_loading_pu": float(st["max_loading_pu"]),
            "voltage_stability_idx": float(st["voltage_stability_idx"]),
            "inertia_mw_s": float(st["inertia_mw_s"]),
            "inverter_tripped_mw": float(st["inverter_tripped_mw"]),
            "frequency_proxy_hz": float(st["frequency_proxy_hz"]),
        }

    def _hourly_state_for_visualization(
        self, hour: int, soc_by_name: Dict[str, float]
    ) -> Dict[str, object]:
        line_service = self._line_service_map()
        st = self._simulate_hour_edge_case(hour, soc_by_name, line_service)
        return {
            "hour": int(st["hour"]),
            "p_inj": st["p_inj"],
            "line_flows": st["line_flows"],
            "max_loading_pu": float(st["max_loading_pu"]),
            "curtailment_mw": float(st["curtailment_mw"]),
            "balancer_mw": float(st["balancer_mw"]),
            "unserved_mw": float(st["unserved_mw"]),
            "emergency_shed_mw": float(st["emergency_shed_mw"]),
            "fault_active": float(st["fault_active"]),
            "fault_type": st["fault_type"],
            "cascaded_lines": st["cascaded_lines"],
            "blackout": float(st["blackout"]),
            "frequency_proxy_hz": float(st["frequency_proxy_hz"]),
            "rocof_proxy_hz_s": float(st["rocof_proxy_hz_s"]),
            "inverter_tripped_mw": float(st["inverter_tripped_mw"]),
            "total_load_mw": float(st["load_mw"]),
            "total_renew_mw": float(st["renewable_mw"]),
        }

    def build_daily_visualization_states(self) -> List[Dict[str, object]]:
        soc_by_name = {s.name: s.soc_mwh for s in self.storage}
        states: List[Dict[str, object]] = []
        for h in range(self.cfg.hours):
            line_service = self._line_service_map()
            st = self._simulate_hour_edge_case(h, soc_by_name, line_service)
            states.append(
                {
                    "hour": int(st["hour"]),
                    "p_inj": st["p_inj"],
                    "line_flows": st["line_flows"],
                    "max_loading_pu": float(st["max_loading_pu"]),
                    "curtailment_mw": float(st["curtailment_mw"]),
                    "balancer_mw": float(st["balancer_mw"]),
                    "unserved_mw": float(st["unserved_mw"]),
                    "emergency_shed_mw": float(st["emergency_shed_mw"]),
                    "fault_active": float(st["fault_active"]),
                    "fault_type": st["fault_type"],
                    "cascaded_lines": st["cascaded_lines"],
                    "blackout": float(st["blackout"]),
                    "frequency_proxy_hz": float(st["frequency_proxy_hz"]),
                    "rocof_proxy_hz_s": float(st["rocof_proxy_hz_s"]),
                    "inverter_tripped_mw": float(st["inverter_tripped_mw"]),
                    "total_load_mw": float(st["load_mw"]),
                    "total_renew_mw": float(st["renewable_mw"]),
                }
            )
        return states

    def animate_load_flow(
        self,
        interval_ms: int = 700,
        save_path: str | None = None,
        show: bool = True,
    ) -> None:
        try:
            import matplotlib.pyplot as plt
            from matplotlib.animation import FuncAnimation
        except ImportError as exc:
            raise RuntimeError(
                "matplotlib is required for animation. Install with: pip install matplotlib"
            ) from exc

        states = self.build_daily_visualization_states()

        bus_xy = {
            0: (0.10, 0.60),
            1: (0.22, 0.82),
            2: (0.42, 0.86),
            3: (0.62, 0.74),
            4: (0.78, 0.54),
            5: (0.64, 0.28),
            6: (0.38, 0.20),
            7: (0.16, 0.34),
        }

        fig = plt.figure(figsize=(14, 8.5))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1.0], height_ratios=[1.0, 1.0])
        ax = fig.add_subplot(gs[:, 0])
        ax_ts1 = fig.add_subplot(gs[0, 1])
        ax_ts2 = fig.add_subplot(gs[1, 1])
        fig.patch.set_facecolor("#f3f6fb")

        ax.set_xlim(0.0, 0.9)
        ax.set_ylim(0.1, 0.95)
        ax.axis("off")
        ax.set_facecolor("#eaf1fb")
        ax.set_title("Renewable Grid Digital Twin: Dynamic Load-Flow", fontsize=14, fontweight="bold")

        line_artists = {}
        for ln in self.lines:
            x1, y1 = bus_xy[ln.from_bus]
            x2, y2 = bus_xy[ln.to_bus]
            (artist,) = ax.plot([x1, x2], [y1, y2], color="gray", linewidth=2.0, alpha=0.8)
            line_artists[ln.name] = artist
            ax.text((x1 + x2) / 2, (y1 + y2) / 2, ln.name, fontsize=8, color="black")

        x = [bus_xy[i][0] for i in range(self.bus_count)]
        y = [bus_xy[i][1] for i in range(self.bus_count)]
        scat = ax.scatter(x, y, s=340, c=["lightgray"] * self.bus_count, edgecolors="black", zorder=5)

        for b in self.buses:
            bx, by = bus_xy[b.idx]
            ax.text(bx, by + 0.04, b.name, ha="center", fontsize=9, fontweight="bold")

        status_text = ax.text(0.02, 0.13, "", fontsize=10, family="monospace")

        hours = [int(s["hour"]) for s in states]
        load_series = [float(s["total_load_mw"]) for s in states]
        renew_series = [float(s["total_renew_mw"]) for s in states]
        bal_series = [float(s["balancer_mw"]) for s in states]
        curtail_series = [float(s["curtailment_mw"]) for s in states]
        line_series = [float(s["max_loading_pu"]) for s in states]
        unserved_series = [float(s["unserved_mw"]) for s in states]
        freq_series = [float(s["frequency_proxy_hz"]) for s in states]
        inv_trip_series = [float(s["inverter_tripped_mw"]) for s in states]
        shed_series = [float(s["emergency_shed_mw"]) for s in states]
        fault_series = [float(s["fault_active"]) for s in states]

        ax_ts1.set_facecolor("white")
        ax_ts1.plot(hours, load_series, color="#b42318", lw=2.0, label="Load MW")
        ax_ts1.plot(hours, renew_series, color="#027a48", lw=2.0, label="Renewables MW")
        ax_ts1.plot(hours, unserved_series, color="#7f1d1d", lw=1.8, label="Unserved MW")
        marker_ts1 = ax_ts1.axvline(0, color="#111827", lw=1.3, ls="--")
        ax_ts1.grid(alpha=0.28)
        ax_ts1.set_title("Power Balance")
        ax_ts1.set_xlabel("Hour")
        ax_ts1.set_ylabel("MW")
        ax_ts1.legend(loc="upper right", fontsize=8)

        ax_ts2.set_facecolor("white")
        ax_ts2.plot(hours, line_series, color="#1d4ed8", lw=2.0, label="Max line loading (pu)")
        ax_ts2.plot(hours, bal_series, color="#c2410c", lw=1.8, label="Balancer MW")
        ax_ts2.plot(hours, curtail_series, color="#6b7280", lw=1.8, label="Curtailment MW")
        ax_ts2.plot(hours, freq_series, color="#6d28d9", lw=1.6, label="Frequency proxy (Hz)")
        ax_ts2.plot(hours, inv_trip_series, color="#be123c", lw=1.4, label="Inverter trip MW")
        ax_ts2.plot(hours, shed_series, color="#0f766e", lw=1.5, label="Emergency shed MW")
        ax_ts2.plot(hours, fault_series, color="#000000", lw=1.2, ls="--", label="Fault flag (0/1)")
        marker_ts2 = ax_ts2.axvline(0, color="#111827", lw=1.3, ls="--")
        ax_ts2.axhline(1.0, color="#b91c1c", lw=1.0, ls=":")
        ax_ts2.grid(alpha=0.28)
        ax_ts2.set_title("Security Stress")
        ax_ts2.set_xlabel("Hour")
        ax_ts2.set_ylabel("pu / MW")
        ax_ts2.legend(loc="upper right", fontsize=8)

        def line_color(loading: float) -> str:
            if loading < 0.70:
                return "#2a9d8f"
            if loading < 1.0:
                return "#e9c46a"
            if loading < 1.25:
                return "#f4a261"
            return "#e63946"

        def update(frame: int):
            st = states[frame]
            flows = st["line_flows"]
            p_inj = st["p_inj"]

            for ln in self.lines:
                f = flows[ln.name]
                loading = abs(f) / ln.thermal_limit_mw
                line_artists[ln.name].set_color(line_color(loading))
                line_artists[ln.name].set_linewidth(1.5 + 5.0 * min(loading, 1.6))

            node_colors = []
            for inj in p_inj:
                if inj > 40:
                    node_colors.append("#2a9d8f")
                elif inj < -40:
                    node_colors.append("#457b9d")
                else:
                    node_colors.append("#bdbdbd")
            if st["blackout"] > 0.5:
                node_colors = ["#7f1d1d" for _ in node_colors]
            scat.set_color(node_colors)

            status_text.set_text(
                f"Hour={st['hour']:02d} | Load={st['total_load_mw']:.0f} MW | "
                f"Renew={st['total_renew_mw']:.0f} MW | "
                f"MaxLine={st['max_loading_pu']:.2f} pu | "
                f"Bal={st['balancer_mw']:.0f} MW | Curt={st['curtailment_mw']:.0f} MW\n"
                f"Fault={st['fault_type']} | Cascades={st['cascaded_lines'] or 'None'} | "
                f"Unserved={st['unserved_mw']:.0f} MW | Freq={st['frequency_proxy_hz']:.2f} Hz | "
                f"Blackout={bool(st['blackout'])}"
            )

            marker_ts1.set_xdata([st["hour"], st["hour"]])
            marker_ts2.set_xdata([st["hour"], st["hour"]])
            return list(line_artists.values()) + [scat, status_text, marker_ts1, marker_ts2]

        ani = FuncAnimation(
            fig,
            update,
            frames=len(states),
            interval=interval_ms,
            blit=False,
            repeat=True,
        )

        if save_path:
            ani.save(save_path, dpi=130)
        if show:
            plt.show()
        else:
            plt.close(fig)

    def _system_inertia(self, hour: int) -> float:
        total = 0.0
        for p in self.plants:
            total += p.inertia_emulation_mw_s * p.profile[hour]
        return total

    def _voltage_proxy(self, p_inj_mw: List[float]) -> float:
        # Fast approximate voltage-risk indicator from net bus stress.
        spread = mean(abs(p) for p in p_inj_mw) / 300.0
        return clamp(1.0 - 0.12 * spread, 0.78, 1.06)

    def run_24h_simulation(self) -> None:
        self.hourly_results.clear()
        soc_by_name = {s.name: s.soc_mwh for s in self.storage}
        for h in range(self.cfg.hours):
            line_service = self._line_service_map()
            st = self._simulate_hour_edge_case(h, soc_by_name, line_service)
            self.hourly_results.append(
                {
                    "hour": float(st["hour"]),
                    "load_mw": float(st["load_mw"]),
                    "renewable_mw": float(st["renewable_mw"]),
                    "storage_net_mw": float(st["storage_net_mw"]),
                    "demand_response_mw": float(st["demand_response_mw"]),
                    "emergency_shed_mw": float(st["emergency_shed_mw"]),
                    "curtailment_mw": float(st["curtailment_mw"]),
                    "balancer_mw": float(st["balancer_mw"]),
                    "unserved_mw": float(st["unserved_mw"]),
                    "hydrogen_load_mw": float(st["hydrogen_load_mw"]),
                    "fault_active": float(st["fault_active"]),
                    "blackout": float(st["blackout"]),
                    "overloaded_lines": float(st["overloaded_lines"]),
                    "max_loading_pu": float(st["max_loading_pu"]),
                    "voltage_stability_idx": float(st["voltage_stability_idx"]),
                    "inertia_mw_s": float(st["inertia_mw_s"]),
                    "inverter_tripped_mw": float(st["inverter_tripped_mw"]),
                    "frequency_proxy_hz": float(st["frequency_proxy_hz"]),
                }
            )

    def contingency_scan(self, target_hour: int = 18) -> List[Dict[str, float]]:
        base = self._hourly_power_balance(target_hour)

        load = [b.load_mw_profile[target_hour] for b in self.buses]
        plant_mw = self._plant_output(target_hour)

        gen_per_bus = [0.0 for _ in range(self.bus_count)]
        for p in self.plants:
            gen_per_bus[p.bus] += plant_mw[p.name]

        p_inj_base = [gen_per_bus[i] - load[i] for i in range(self.bus_count)]

        findings = []
        for idx, tripped in enumerate(self.lines):
            case_lines = [
                Line(
                    name=l.name,
                    from_bus=l.from_bus,
                    to_bus=l.to_bus,
                    x_pu=l.x_pu,
                    thermal_limit_mw=l.thermal_limit_mw,
                    in_service=(l.in_service and i != idx),
                )
                for i, l in enumerate(self.lines)
            ]
            try:
                _, flows = self.dc_load_flow(p_inj_base, case_lines)
                max_lp = max(abs(flows[l.name]) / l.thermal_limit_mw for l in self.lines)
                findings.append(
                    {
                        "outage": tripped.name,
                        "islanded": 0.0,
                        "max_loading_pu": max_lp,
                        "severity": max(0.0, max_lp - 1.0),
                    }
                )
            except ValueError:
                findings.append(
                    {
                        "outage": tripped.name,
                        "islanded": 1.0,
                        "max_loading_pu": 9.9,
                        "severity": 9.9,
                    }
                )

        findings.sort(key=lambda x: x["severity"], reverse=True)
        findings.insert(
            0,
            {
                "outage": "BASE",
                "islanded": 0.0,
                "max_loading_pu": base["max_loading_pu"],
                "severity": max(0.0, base["max_loading_pu"] - 1.0),
            },
        )
        return findings

    def transient_frequency_response(
        self,
        disturbance_mw: float,
        hour: int = 19,
        damping_pu: float = 1.4,
    ) -> Dict[str, float]:
        # Simplified center-of-inertia frequency model:
        # df/dt = (DeltaP - D*df) / (2H)
        # where H is represented from aggregated inertia MW*s.
        f0 = self.cfg.nominal_frequency_hz
        dt = self.cfg.dt_s
        steps = int(self.cfg.transient_window_s / dt)

        inertia = max(self._system_inertia(hour), 80.0)
        synthetic_fast_freq = 0.15 * disturbance_mw

        freq = f0
        f_min = f0
        rocof_max = 0.0

        for _ in range(steps):
            effective_delta_p = disturbance_mw - synthetic_fast_freq
            dfdot = -(effective_delta_p + damping_pu * (f0 - freq) * 100) / (2.0 * inertia)
            freq += dfdot * dt
            f_min = min(f_min, freq)
            rocof_max = max(rocof_max, abs(dfdot))

        return {
            "disturbance_mw": disturbance_mw,
            "inertia_mw_s": inertia,
            "frequency_nadir_hz": f_min,
            "steady_frequency_hz": freq,
            "max_rocof_hz_per_s": rocof_max,
            "secure": float(f_min >= 49.0 and rocof_max < 1.0),
        }

    def monte_carlo_adequacy(self, samples: int = 200) -> Dict[str, float]:
        lolp_count = 0
        eens_total = 0.0

        for _ in range(samples):
            h = self.rng.randrange(0, self.cfg.hours)
            weather_scale = clamp(self.rng.gauss(1.0, 0.2), 0.5, 1.3)
            forced_outage = 1 if self.rng.random() < 0.07 else 0

            generation = 0.0
            for p in self.plants:
                p_mw = p.p_max_mw * p.profile[h] * weather_scale
                if forced_outage and p.plant_type in {"wind", "solar"} and self.rng.random() < 0.2:
                    p_mw *= 0.3
                generation += p_mw

            load = sum(b.load_mw_profile[h] for b in self.buses) * clamp(self.rng.gauss(1.0, 0.05), 0.9, 1.15)
            deficit = max(0.0, load - generation)
            if deficit > 0:
                lolp_count += 1
            eens_total += deficit

        return {
            "samples": float(samples),
            "LOLP": lolp_count / samples,
            "EENS_MWh_per_sample": eens_total / samples,
        }

    def edge_stability_scan(self, top_n: int = 5) -> List[Dict[str, float]]:
        local_soc = {s.name: s.soc_mwh for s in self.storage}
        local_lines = self._line_service_map()
        records: List[Dict[str, float]] = []

        for hour in range(self.cfg.hours):
            st = self._simulate_hour_edge_case(hour, local_soc, local_lines)
            severity = (
                2.2 * float(st["max_loading_pu"])
                + 0.03 * float(st["unserved_mw"])
                + 9.0 * max(0.0, 49.0 - float(st["frequency_proxy_hz"]))
                + 3.0 * float(st["blackout"])
                + 1.7 * float(st["fault_active"])
            )
            records.append(
                {
                    "hour": float(hour),
                    "severity": severity,
                    "max_loading_pu": float(st["max_loading_pu"]),
                    "unserved_mw": float(st["unserved_mw"]),
                    "frequency_proxy_hz": float(st["frequency_proxy_hz"]),
                    "blackout": float(st["blackout"]),
                    "fault_active": float(st["fault_active"]),
                }
            )
        records.sort(key=lambda r: r["severity"], reverse=True)
        return records[:top_n]

    def summary(self) -> str:
        if not self.hourly_results:
            return "No simulation data. Run run_24h_simulation() first."

        avg_renew = mean(r["renewable_mw"] for r in self.hourly_results)
        avg_load = mean(r["load_mw"] for r in self.hourly_results)
        avg_curtail = mean(r["curtailment_mw"] for r in self.hourly_results)
        avg_unserved = mean(r["unserved_mw"] for r in self.hourly_results)
        avg_freq = mean(r["frequency_proxy_hz"] for r in self.hourly_results)
        total_blackouts = int(sum(1 for r in self.hourly_results if r["blackout"] > 0.5))
        total_fault_hours = int(sum(1 for r in self.hourly_results if r["fault_active"] > 0.5))
        total_inverter_trip = sum(r["inverter_tripped_mw"] for r in self.hourly_results)
        peak_stress = max(self.hourly_results, key=lambda r: r["max_loading_pu"])
        min_voltage = min(r["voltage_stability_idx"] for r in self.hourly_results)

        lines = [
            "=== Renewable-Dominant Power System Study ===",
            f"Buses: {self.bus_count}, Lines: {len(self.lines)}, Plants: {len(self.plants)}, BESS: {len(self.storage)}",
            f"Average Load: {avg_load:.1f} MW",
            f"Average Renewable Output: {avg_renew:.1f} MW",
            f"Average Curtailment: {avg_curtail:.1f} MW",
            f"Average Unserved Energy: {avg_unserved:.1f} MW",
            f"Average Frequency Proxy: {avg_freq:.2f} Hz",
            f"Peak Network Loading: {peak_stress['max_loading_pu']:.2f} pu at hour {int(peak_stress['hour'])}",
            f"Minimum Voltage Stability Proxy: {min_voltage:.3f} pu",
            f"Fault Hours: {total_fault_hours}, Blackout Hours: {total_blackouts}",
            f"Total Inverter Trip Energy-Loss Proxy: {total_inverter_trip:.1f} MW",
        ]
        return "\n".join(lines)


def run_demo(
    show_animation: bool = False,
    save_animation_path: str | None = None,
    interval_ms: int = 750,
) -> None:
    model = RenewableGridModel(SimulationConfig())
    model.run_24h_simulation()

    print(model.summary())

    n1 = model.contingency_scan(target_hour=18)
    print("\nTop 5 N-1 contingencies (higher severity is worse):")
    for row in n1[:5]:
        print(
            f"  {row['outage']:<5} | islanded={int(row['islanded'])} | "
            f"max_loading={row['max_loading_pu']:.2f} pu | severity={row['severity']:.2f}"
        )

    freq = model.transient_frequency_response(disturbance_mw=220.0, hour=19)
    print("\nTransient frequency stability:")
    print(
        f"  Nadir={freq['frequency_nadir_hz']:.3f} Hz, "
        f"RoCoF={freq['max_rocof_hz_per_s']:.3f} Hz/s, "
        f"Secure={bool(freq['secure'])}"
    )

    adequacy = model.monte_carlo_adequacy(samples=300)
    print("\nMonte Carlo adequacy:")
    print(
        f"  LOLP={adequacy['LOLP']:.3f}, "
        f"EENS={adequacy['EENS_MWh_per_sample']:.2f} MWh/sample"
    )

    edge = model.edge_stability_scan(top_n=5)
    print("\nEdge stability critical windows:")
    for row in edge:
        print(
            f"  h={int(row['hour']):02d} | severity={row['severity']:.2f} | "
            f"maxLine={row['max_loading_pu']:.2f} pu | unserved={row['unserved_mw']:.1f} MW | "
            f"f={row['frequency_proxy_hz']:.2f} Hz | blackout={bool(row['blackout'])}"
        )

    if show_animation or save_animation_path:
        model.animate_load_flow(
            interval_ms=interval_ms,
            save_path=save_animation_path,
            show=show_animation,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Renewable grid simulation with animated load-flow dashboard.")
    parser.add_argument("--show", action="store_true", help="Open interactive animation window.")
    parser.add_argument("--save-gif", type=str, default=None, help="Save animation GIF path.")
    parser.add_argument("--interval-ms", type=int, default=750, help="Animation frame interval in milliseconds.")
    parser.add_argument(
        "--auto-fallback-gif",
        type=str,
        default="renewable_load_flow.gif",
        help="Auto-save GIF path when backend is non-interactive and --show is selected.",
    )
    args = parser.parse_args()

    show_requested = args.show or (args.save_gif is None)
    save_path = args.save_gif

    try:
        import matplotlib.pyplot as plt

        backend_name = str(plt.get_backend()).lower().strip()
        non_interactive_backends = {"agg", "pdf", "svg", "ps", "template", "cairo"}
        backend_non_interactive = backend_name in non_interactive_backends
    except Exception:
        backend_non_interactive = False

    if show_requested and backend_non_interactive and save_path is None:
        save_path = args.auto_fallback_gif
        show_requested = False
        print(
            f"Non-interactive backend detected. "
            f"Saving dashboard animation to '{save_path}' instead of opening a window."
        )

    try:
        run_demo(
            show_animation=show_requested,
            save_animation_path=save_path,
            interval_ms=args.interval_ms,
        )
    except Exception as exc:
        fallback = save_path or "renewable_load_flow.gif"
        print(f"\nAnimation display failed: {exc}")
        print(f"Trying to save animation to '{fallback}' instead...")
        run_demo(show_animation=False, save_animation_path=fallback, interval_ms=args.interval_ms)
        print(f"Saved animation to {fallback}")

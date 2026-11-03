from pathlib import Path

html = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Simplex LP Solver</title>
<style>
  :root {
    --bg:#0b1020; --panel:#11182b; --panel2:#18223a; --text:#e8edf7;
    --muted:#9aa7bd; --accent:#67b7ff; --good:#49d17d; --bad:#ff6b6b;
    --border:#2a3857; --input:#0d1426;
  }
  *{box-sizing:border-box}
  body{margin:0;font-family:Inter,system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--text)}
  .wrap{max-width:1500px;margin:auto;padding:24px}
  h1{margin:0 0 6px;font-size:30px}
  h2{font-size:18px;margin:0 0 14px}
  .sub{color:var(--muted);margin-bottom:22px}
  .grid{display:grid;grid-template-columns:360px 1fr;gap:18px}
  .card{background:var(--panel);border:1px solid var(--border);border-radius:14px;padding:18px}
  .controls{display:grid;grid-template-columns:1fr 1fr;gap:10px}
  label{font-size:13px;color:var(--muted);display:block;margin-bottom:5px}
  input,select,button{font:inherit}
  input,select{width:100%;background:var(--input);color:var(--text);border:1px solid var(--border);border-radius:8px;padding:9px}
  input:focus,select:focus{outline:2px solid #2d76b7;border-color:#67b7ff}
  button{border:0;border-radius:9px;padding:10px 13px;cursor:pointer;background:var(--panel2);color:var(--text);border:1px solid var(--border)}
  button:hover{filter:brightness(1.15)}
  button.primary{background:#1769aa;border-color:#2788d4}
  button.danger{background:#59212a;border-color:#87303d}
  .btnrow{display:flex;gap:8px;flex-wrap:wrap;margin-top:12px}
  .objective{overflow-x:auto;margin-top:14px}
  .objrow,.constraint{display:flex;gap:7px;align-items:center;margin-bottom:8px;min-width:max-content}
  .objrow input{width:82px}
  .constraint input{width:82px}
  .constraint .sign{width:72px}
  .constraint .rhs{width:95px}
  .varname{color:var(--accent);font-weight:600;min-width:25px}
  .sectionhead{display:flex;justify-content:space-between;align-items:center;gap:10px;margin-bottom:12px}
  .hint{font-size:12px;color:var(--muted);line-height:1.5}
  .status{padding:12px;border-radius:10px;margin-top:14px;background:var(--panel2);border-left:4px solid var(--accent)}
  .status.good{border-left-color:var(--good)}
  .status.bad{border-left-color:var(--bad)}
  .results{display:none}
  .resultbox{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:9px;margin-bottom:16px}
  .metric{background:var(--panel2);padding:12px;border-radius:10px}
  .metric small{display:block;color:var(--muted);margin-bottom:4px}
  .metric b{font-size:18px}
  .tablewrap{overflow:auto;border:1px solid var(--border);border-radius:10px;margin-top:10px}
  table{border-collapse:collapse;width:100%;min-width:650px}
  th,td{padding:7px 9px;border-bottom:1px solid var(--border);text-align:right;white-space:nowrap;font-size:13px}
  th{background:#1a2742;position:sticky;top:0;z-index:1}
  th:first-child,td:first-child{text-align:left}
  .step{margin-top:18px}
  .step h3{font-size:14px;color:var(--accent);margin:0 0 7px}
  .note{background:#0e172a;border:1px solid var(--border);padding:10px;border-radius:8px;color:var(--muted);font-size:12px}
  .check{display:flex;align-items:center;gap:8px;margin-top:10px;color:var(--muted);font-size:13px}
  .check input{width:auto}
  @media(max-width:900px){.grid{grid-template-columns:1fr}.wrap{padding:14px}}
</style>
</head>
<body>
<div class="wrap">
  <h1>Simplex Linear Programming Solver</h1>
  <div class="sub">Arbitrary number of variables and constraints • ≤, ≥ and = constraints • maximization/minimization • two-phase simplex</div>

  <div class="grid">
    <div class="card">
      <div class="sectionhead"><h2>Problem setup</h2></div>

      <div class="controls">
        <div>
          <label>Number of variables</label>
          <input id="nVars" type="number" min="1" max="100" value="2">
        </div>
        <div>
          <label>Number of constraints</label>
          <input id="nCons" type="number" min="1" max="100" value="2">
        </div>
      </div>

      <div class="btnrow">
        <button onclick="buildProblem()">Apply dimensions</button>
        <button onclick="loadExample()">Load example</button>
      </div>

      <div class="objective">
        <h2 style="margin-top:18px">Objective function</h2>
        <div class="controls" style="grid-template-columns:120px 1fr">
          <div>
            <label>Direction</label>
            <select id="direction">
              <option value="max">Maximize</option>
              <option value="min">Minimize</option>
            </select>
          </div>
          <div>
            <label>Coefficients</label>
            <div id="objectiveInputs" class="objrow"></div>
          </div>
        </div>
      </div>

      <h2 style="margin-top:18px">Constraints</h2>
      <div id="constraints"></div>

      <label class="check">
        <input id="showSteps" type="checkbox" checked>
        Show simplex iterations
      </label>

      <div class="btnrow">
        <button class="primary" onclick="solve()">Solve</button>
        <button onclick="clearInputs()">Clear values</button>
      </div>

      <div id="status" class="status">Enter your LP and press <b>Solve</b>.</div>
      <div class="hint" style="margin-top:12px">
        Assumption: all decision variables are non-negative (xᵢ ≥ 0).
        For a free/unrestricted variable, substitute x = x⁺ − x⁻.
      </div>
    </div>

    <div class="card">
      <div class="sectionhead"><h2>Solution</h2></div>
      <div id="results" class="results">
        <div id="metrics" class="resultbox"></div>
        <div id="solutionTable"></div>
        <div id="steps"></div>
      </div>
      <div id="empty" class="note">
        Your optimal solution and tableau iterations will appear here.
      </div>
    </div>
  </div>
</div>

<script>
const EPS = 1e-9;
const MAX_ITER = 5000;

function fmt(x) {
  if (!Number.isFinite(x)) return "∞";
  if (Math.abs(x) < 1e-8) x = 0;
  return Number(x).toFixed(6).replace(/\.?0+$/,'');
}
function val(id) { return Number(document.getElementById(id).value); }

function buildProblem() {
  const n = Math.max(1, Math.min(100, parseInt(document.getElementById('nVars').value)||1));
  const m = Math.max(1, Math.min(100, parseInt(document.getElementById('nCons').value)||1));
  document.getElementById('nVars').value=n;
  document.getElementById('nCons').value=m;

  const oi = document.getElementById('objectiveInputs');
  oi.innerHTML='';
  for(let j=0;j<n;j++){
    const inp=document.createElement('input');
    inp.type='number'; inp.step='any'; inp.value=j===0?'3':j===1?'2':'0';
    inp.id='c_'+j; inp.title='Coefficient of x'+(j+1);
    oi.appendChild(inp);
    const span=document.createElement('span');
    span.textContent=' x'+(j+1)+(j<n-1?' + ':'');
    oi.appendChild(span);
  }

  const box=document.getElementById('constraints');
  box.innerHTML='';
  for(let i=0;i<m;i++){
    const row=document.createElement('div'); row.className='constraint'; row.id='con_'+i;
    for(let j=0;j<n;j++){
      const inp=document.createElement('input');
      inp.type='number'; inp.step='any'; inp.value=(i===0&&j===0)?1:(i===0&&j===1)?1:(i===1&&j===0)?2:(i===1&&j===1)?1:0;
      inp.id=`a_${i}_${j}`; row.appendChild(inp);
      const sp=document.createElement('span'); sp.textContent=' x'+(j+1)+(j<n-1?' + ':''); row.appendChild(sp);
    }
    const sel=document.createElement('select'); sel.className='sign'; sel.id='sign_'+i;
    sel.innerHTML='<option value="<=">≤</option><option value=">=">≥</option><option value="=">=</option>';
    row.appendChild(sel);
    const rhs=document.createElement('input'); rhs.type='number'; rhs.step='any'; rhs.className='rhs'; rhs.id='b_'+i;
    rhs.value=i===0?4:i===1?5:0; row.appendChild(rhs);
    box.appendChild(row);
  }
}
function clearInputs(){
  document.querySelectorAll('input[type=number]').forEach(x=>x.value='');
  document.getElementById('results').style.display='none';
  document.getElementById('empty').style.display='block';
  setStatus('Inputs cleared.','');
}
function loadExample(){
  document.getElementById('nVars').value=2;
  document.getElementById('nCons').value=2;
  buildProblem();
  document.getElementById('direction').value='max';
  document.getElementById('c_0').value=3; document.getElementById('c_1').value=2;
  document.getElementById('a_0_0').value=1; document.getElementById('a_0_1').value=1; document.getElementById('sign_0').value='<='; document.getElementById('b_0').value=4;
  document.getElementById('a_1_0').value=2; document.getElementById('a_1_1').value=1; document.getElementById('sign_1').value='<='; document.getElementById('b_1').value=5;
}
function setStatus(msg, cls){
  const s=document.getElementById('status'); s.innerHTML=msg; s.className='status '+cls;
}

function readProblem(){
  const n=Number(document.getElementById('nVars').value);
  const m=Number(document.getElementById('nCons').value);
  if(!Number.isInteger(n)||n<1||n>100||!Number.isInteger(m)||m<1||m>100) throw new Error('Use 1–100 variables and 1–100 constraints.');
  const c=[];
  for(let j=0;j<n;j++){
    const x=Number(document.getElementById('c_'+j).value);
    if(!Number.isFinite(x)) throw new Error('Enter every objective coefficient.');
    c.push(x);
  }
  const A=[], b=[], signs=[];
  for(let i=0;i<m;i++){
    const row=[];
    for(let j=0;j<n;j++){
      const x=Number(document.getElementById(`a_${i}_${j}`).value);
      if(!Number.isFinite(x)) throw new Error(`Missing coefficient in constraint ${i+1}.`);
      row.push(x);
    }
    const rhs=Number(document.getElementById('b_'+i).value);
    if(!Number.isFinite(rhs)) throw new Error(`Missing RHS in constraint ${i+1}.`);
    A.push(row); b.push(rhs); signs.push(document.getElementById('sign_'+i).value);
  }
  return {n,m,c,A,b,signs,direction:document.getElementById('direction').value};
}

function standardize(p){
  // Build rows with slack/surplus/artificial variables.
  let names=p.c.map((_,j)=>'x'+(j+1));
  let rows=[], rhs=[], basis=[], artificial=new Set();
  for(let i=0;i<p.m;i++){
    let row=p.A[i].slice(), r=p.b[i], sign=p.signs[i];
    // Make RHS nonnegative. Reversing inequality also reverses its sense.
    if(r < -EPS){
      row=row.map(v=>-v); r=-r;
      if(sign==='<=' ) sign='>=';
      else if(sign==='>=') sign='<=';
    }
    while(row.length<names.length) row.push(0);
    if(sign==='<='){
      const idx=names.length; names.push('s'+(i+1));
      row.push(1); basis.push(idx);
    } else if(sign==='>='){
      const sur=names.length; names.push('u'+(i+1)); row.push(-1);
      const art=names.length; names.push('a'+(i+1)); row.push(1);
      basis.push(art); artificial.add(art);
    } else {
      const art=names.length; names.push('a'+(i+1)); row.push(1);
      basis.push(art); artificial.add(art);
    }
    rows.push(row); rhs.push(r);
  }
  // Pad all rows to equal length.
  const N=names.length;
  rows=rows.map(r=>r.concat(Array(N-r.length).fill(0)));
  return {names,rows,rhs,basis,artificial};
}

function makeTableau(std, costs){
  const m=std.rows.length, N=std.names.length;
  const T=std.rows.map((r,i)=>r.slice().concat([std.rhs[i]]));
  const z=Array(N+1).fill(0);
  for(let j=0;j<N;j++) z[j]=-costs[j];
  // Canonicalize objective row: eliminate basic variables.
  for(let i=0;i<m;i++){
    const b=std.basis[i], cb=costs[b]||0;
    if(Math.abs(cb)>EPS){
      for(let j=0;j<=N;j++) z[j]+=cb*T[i][j];
    }
  }
  T.push(z);
  return T;
}
function pivot(T,row,col){
  const m=T.length, n=T[0].length;
  const pv=T[row][col];
  for(let j=0;j<n;j++) T[row][j]/=pv;
  for(let i=0;i<m;i++){
    if(i===row) continue;
    const f=T[i][col];
    if(Math.abs(f)>EPS){
      for(let j=0;j<n;j++) T[i][j]-=f*T[row][j];
    }
  }
}
function simplex(T,basis,nVars,allowedCols,steps,phaseLabel,names){
  const m=T.length-1, N=T[0].length-1;
  let iter=0;
  while(iter<MAX_ITER){
    let enter=-1, most=EPS;
    for(const j of allowedCols){
      if(T[m][j] < -most){ most=-T[m][j]; enter=j; }
    }
    if(enter<0) return {status:'optimal',iterations:iter};
    let leave=-1, best=Infinity;
    for(let i=0;i<m;i++){
      const a=T[i][enter];
      if(a>EPS){
        const ratio=T[i][N]/a;
        if(ratio<best-EPS || (Math.abs(ratio-best)<=EPS && (leave<0 || basis[i]>basis[leave]))){
          best=ratio; leave=i;
        }
      }
    }
    if(leave<0) return {status:'unbounded',iterations:iter};
    steps.push({label:`${phaseLabel} — iteration ${iter+1}`,T:T.map(r=>r.slice()),basis:basis.slice(),names:names.slice(),enter,leave});
    pivot(T,leave,enter);
    basis[leave]=enter;
    iter++;
  }
  return {status:'iteration_limit',iterations:iter};
}

function removeArtificial(T,basis,names,artificial){
  let m=T.length-1, N=T[0].length-1;
  // Pivot artificial basic variables out where possible; otherwise remove redundant zero rows.
  for(let i=0;i<m;i++){
    if(!artificial.has(basis[i])) continue;
    let col=-1;
    for(let j=0;j<N;j++){
      if(!artificial.has(j) && Math.abs(T[i][j])>EPS){ col=j; break; }
    }
    if(col>=0){
      pivot(T,i,col); basis[i]=col;
    } else {
      // Redundant row if RHS is zero.
      if(Math.abs(T[i][N])<EPS){
        T.splice(i,1); basis.splice(i,1); m--; i--;
      } else {
        throw new Error('Infeasible: an artificial basic variable remains positive.');
      }
    }
  }
  // Delete artificial columns.
  const keep=[];
  for(let j=0;j<N;j++) if(!artificial.has(j)) keep.push(j);
  const map=new Map(keep.map((old,k)=>[old,k]));
  const newT=T.map(r=>{
    const nr=keep.map(j=>r[j]);
    nr.push(r[N]);
    return nr;
  });
  const newBasis=basis.map(b=>map.get(b));
  const newNames=keep.map(j=>names[j]);
  return {T:newT,basis:newBasis,names:newNames};
}

function solveLP(p){
  const std=standardize(p);
  let cMax=p.c.slice();
  if(p.direction==='min') cMax=cMax.map(x=>-x);
  const costs=Array(std.names.length).fill(0);
  for(let j=0;j<p.n;j++) costs[j]=cMax[j];

  const steps=[];
  let T,basis,names;
  // Phase I if artificial variables exist.
  if(std.artificial.size){
    const phase1Costs=Array(std.names.length).fill(0);
    for(const a of std.artificial) phase1Costs[a]=-1;
    T=makeTableau(std,phase1Costs);
    basis=std.basis.slice(); names=std.names.slice();
    const artCols=[...Array(names.length).keys()];
    const r1=simplex(T,basis,p.n,artCols,steps,'Phase I',names);
    if(r1.status==='unbounded') return {status:'infeasible',message:'Phase I became unbounded; the model is infeasible or numerically problematic.'};
    const phase1Value=T[T.length-1][T[0].length-1];
    if(phase1Value < -1e-7) return {status:'infeasible',message:'No feasible solution exists (Phase I optimum is negative).'};
    const cleaned=removeArtificial(T,basis,names,std.artificial);
    T=cleaned.T; basis=cleaned.basis; names=cleaned.names;
  } else {
    T=std.rows.map((r,i)=>r.slice().concat([std.rhs[i]]));
    basis=std.basis.slice(); names=std.names.slice();
    T.push(Array(T[0].length).fill(0));
  }

  // Phase II objective after artificial columns have been removed.
  const N=names.length;
  const c2=Array(N).fill(0);
  for(let j=0;j<N;j++){
    const originalIndex = names[j].startsWith('x') ? Number(names[j].slice(1))-1 : -1;
    if(originalIndex>=0) c2[j]=cMax[originalIndex];
  }
  const z=Array(N+1).fill(0);
  for(let j=0;j<N;j++) z[j]=-c2[j];
  for(let i=0;i<basis.length;i++){
    const cb=c2[basis[i]]||0;
    if(Math.abs(cb)>EPS) for(let j=0;j<=N;j++) z[j]+=cb*T[i][j];
  }
  T[T.length-1]=z;

  const allowed=[...Array(N).keys()];
  const r2=simplex(T,basis,p.n,allowed,steps,'Phase II',names);
  if(r2.status==='unbounded') return {status:'unbounded',message:'The objective can increase without bound.'};
  if(r2.status==='iteration_limit') return {status:'iteration_limit',message:'Iteration limit reached.'};

  const x=Array(p.n).fill(0);
  for(let i=0;i<basis.length;i++){
    const bvar=basis[i];
    if(names[bvar] && names[bvar].startsWith('x')){
      const k=Number(names[bvar].slice(1))-1;
      if(k>=0&&k<p.n) x[k]=Math.abs(T[i][N])<EPS?0:T[i][N];
    }
  }
  const originalZ=p.c.reduce((s,c,j)=>s+c*x[j],0);
  return {status:'optimal',x,objective:originalZ,T,basis,names,steps};
}

function renderTable(step){
  const T=step.T, names=step.names, basis=step.basis, N=names.length, m=T.length-1;
  let h='<div class="tablewrap"><table><thead><tr><th>Basis</th>';
  for(const n of names) h+=`<th>${n}</th>`;
  h+='<th>RHS</th></tr></thead><tbody>';
  for(let i=0;i<m;i++){
    h+=`<tr><td>${names[basis[i]]||'—'}</td>`;
    for(let j=0;j<=N;j++) h+=`<td>${fmt(T[i][j])}</td>`;
    h+='</tr>';
  }
  h+='<tr><td><b>Z</b></td>';
  for(let j=0;j<=N;j++) h+=`<td><b>${fmt(T[m][j])}</b></td>`;
  h+='</tr></tbody></table></div>';
  return h;
}

function solve(){
  try{
    const p=readProblem();
    const result=solveLP(p);
    const res=document.getElementById('results'), empty=document.getElementById('empty');
    res.style.display='block'; empty.style.display='none';

    if(result.status!=='optimal'){
      document.getElementById('metrics').innerHTML='';
      document.getElementById('solutionTable').innerHTML=`<div class="status bad"><b>${result.status.toUpperCase()}</b><br>${result.message}</div>`;
      document.getElementById('steps').innerHTML='';
      setStatus(result.message,'bad');
      return;
    }

    setStatus('Optimal solution found.','good');
    let metrics=`<div class="metric"><small>Objective</small><b>${fmt(result.objective)}</b></div>`;
    metrics+=`<div class="metric"><small>Variables</small><b>${p.n}</b></div>`;
    metrics+=`<div class="metric"><small>Constraints</small><b>${p.m}</b></div>`;
    metrics+=`<div class="metric"><small>Iterations</small><b>${result.steps.length}</b></div>`;
    document.getElementById('metrics').innerHTML=metrics;

    let sh='<h2>Decision variables</h2><div class="tablewrap"><table><thead><tr><th>Variable</th><th>Value</th></tr></thead><tbody>';
    result.x.forEach((v,j)=>sh+=`<tr><td>x${j+1}</td><td>${fmt(v)}</td></tr>`);
    sh+='</tbody></table></div>';
    document.getElementById('solutionTable').innerHTML=sh;

    const stepsBox=document.getElementById('steps');
    if(!document.getElementById('showSteps').checked){
      stepsBox.innerHTML='<div class="note" style="margin-top:15px">Iterations hidden. Enable “Show simplex iterations” and solve again.</div>';
    } else {
      let out='<div class="step"><h2>Simplex iterations</h2>';
      if(result.steps.length===0) out+='<div class="note">The initial basic feasible solution was already optimal.</div>';
      result.steps.forEach((s,i)=>{
        out+=`<div class="step"><h3>${i+1}. ${s.label}</h3>${renderTable(s)}</div>`;
      });
      out+='</div>';
      stepsBox.innerHTML=out;
    }
  }catch(e){
    setStatus(e.message,'bad');
  }
}
buildProblem();
</script>
</body>
</html>'''

path = Path("/mnt/data/simplex_lp_solver.html")
path.write_text(html, encoding="utf-8")
print(path)

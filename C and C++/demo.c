#include<stdio.h>

struct point
{
    int x;
    int y;
};

struct point *createpoint(int x, int y)
{
    struct point newpoint = {x, y};
    return &newpoint;

};

void main(){
    struct point *p = createpoint(10, 20);
    printf("%d %d", p->x, p->y);
}

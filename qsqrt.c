#include <stdio.h>
float OE_Q_sqrt( float number )
{
    long i;
    float x2, y;
    int fix = 0;
    const float threehalfs = 1.5F;
    x2 = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;                       // evil floating point bit level hacking
    i  = 0x1fc00000 + fix + ( i >> 1 );               // what the fuck? 
    y  = * ( float * ) &i;
    printf("该数字的粗略平方根为： %f\n", y);
    y  = y - (y*y-number)/(2*y);   // 1st iteration
//  y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
    return y;
}
int main() {
    float num;
    printf("请输入一个浮点数：\n");
    scanf("%f", &num); 
    float ret = OE_Q_sqrt(num);
    printf("经过一次牛顿迭代法该数字平方根为： %f\n", y);
    return 0;
}
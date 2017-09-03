#include <stdio.h>

#define ACNT 1000000
#define BCNT 100
#define BCNT2 200
#define BCNT3 400

double cascsum(double p[], int p_len){
    double e = 0;
    double s = p[0];
    double x, z, y;
    for(int j=1; j<p_len; j++){
        x = p[j] + s;
        z = x - p[j];
        y = (p[j] - (x - z)) + (s - z);
        e += y;
        s = x;
    }
    return s + e;
}

double sum(double p[], int p_len){
    double s = p[0];
    for(int j=1; j<p_len; j++){
        s += p[j];
    }
    return s;
}

int main(){
    double num1[ACNT];
    double num2[BCNT3];
    for(int i=0; i<ACNT; i++){
        num1[i] = 1.0/ACNT;
    }
    num2[0] = 1;
    for(int i=1; i<BCNT3; i++){
        num2[i] = num2[i-1]*.99;
    }
    
    //Second series to BCNT terms
    double exact1 = 1.0;
    for(int i=1; i<BCNT+1; i++){
        exact1 *= .99;
    }
    exact1 = (1-exact1)/(1-.99);
    
    //Second series to BCNT2 terms
    double exact2 = 1.0;
    for(int i=1; i<BCNT2+1; i++){
        exact2 *= .99;
    }
    exact2 = (1-exact2)/(1-.99);
    
    //Second series to BCNT3 terms
    double exact3 = 1.0;
    for(int i=1; i<BCNT3+1; i++){
        exact3 *= .99;
    }
    exact3 = (1-exact3)/(1-.99);
    
    printf("For the first series:\n");
    printf("The real sum should be 1.\n");
    printf("Sum     %.16f\n", sum(num1,ACNT));
    printf("Cascsum %.16f\n", cascsum(num1,ACNT));
    printf("Sum error:    %.16e\n", 1.0-sum(num1,ACNT));
    printf("Cascsum error %.16e\n\n", 1.0-cascsum(num1,ACNT));
    
    printf("For the second series:\n");
    printf("The real sum should be\n        %.20f\n", exact1);
    printf("Sum     %.20f\n", sum(num2,BCNT));
    printf("Cascsum %.20f\n", cascsum(num2,BCNT));
    printf("Sum error:    %.16e\n", exact1-sum(num2,BCNT));
    printf("Cascsum error %.16e\n\n", exact1-cascsum(num2,BCNT));
    
    printf("For the second series taken to %d terms:\n", BCNT2);
    printf("The real sum should be\n        %.20f\n", exact2);
    printf("Sum     %.20f\n", sum(num2,BCNT2));
    printf("Cascsum %.20f\n", cascsum(num2,BCNT2));
    printf("Sum error:    %.16e\n", exact2-sum(num2,BCNT2));
    printf("Cascsum error %.16e\n\n", exact2-cascsum(num2,BCNT2));
    
    printf("For the second series taken to %d terms:\n", BCNT3);
    printf("The real sum should be\n        %.20f\n", exact3);
    printf("Sum     %.20f\n", sum(num2,BCNT3));
    printf("Cascsum %.20f\n", cascsum(num2,BCNT3));
    printf("Sum error:    %.16e\n", exact3-sum(num2,BCNT3));
    printf("Cascsum error %.16e\n\n", exact3-cascsum(num2,BCNT3));
}

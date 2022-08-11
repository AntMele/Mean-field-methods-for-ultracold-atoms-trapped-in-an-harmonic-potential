#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int n = 1000;
double h = 0.01;
double Na = -0.574;
double alpha = 0.1;
double * R; // to be used in the potential


double k_square(double E, int i){
    return 2*E - pow(h*i,2) - 2*Na*R[i]/pow((h*i),2);
}

double numerov(double E){
    double k_square0 = 0;
    double k_square1 = k_square(E, 1);
    double k_square2 = k_square(E, 2);
    double y0 = 0;
    double y1 = h;
    double y2 = y1 * (2-5*pow(h,2)*k_square1/6) / (1+pow(h,2)*k_square2/12);
    for(int i = 3; i <= n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square(E, i);
        y0 = y1;
        y1 = y2;
        y2 = (y1 * (2-5*pow(h,2)*k_square1/6) - y0 * (1+pow(h,2)*k_square0/12))/ (1+pow(h,2)*k_square2/12);
    }
    return y2;
}

double numerovSave(double E, double * y){
    double k_square0 = 0;
    double k_square1 = k_square(E, 1);
    double k_square2 = k_square(E, 2);
    y[0] = 0;
    y[1] = h;
    y[2] = y[1] * (2-5*pow(h,2)*k_square1/6) / (1+pow(h,2)*k_square2/12);
    for(int i = 2; i < n; i++){
        k_square0 = k_square1;
        k_square1 = k_square2;
        k_square2 = k_square(E, i+1);
        y[i+1] = (y[i] * (2-5*pow(h,2)*k_square1/6) - y[i-1] * (1+pow(h,2)*k_square0/12))/ (1+pow(h,2)*k_square2/12);
    }

    // I correct the divergence of the wave function at rmax
    double min = fabs(y[1]);
    int index = 1;
    for(int i = 2; i <= n; i++){
        if(fabs(y[i]) < min){
            min = fabs(y[i]);
            index = i;
        }
    }
    for(int i = index+1; i <= n; i++){
        y[i] = min;
    }

    // I calculate the norm with the trapeziudal rule
    double norm = (pow(y[0],2) + pow(y[n],2)) / 2;
    for(int i = 1; i < n; i++){
        norm += pow(y[i], 2);
    }
    norm = sqrt(norm * h);
    for(int i = 0; i <= n; i++){
        y[i] /= norm;
    }

    return y[n];
}

double secantMethod(double xmin, double xmax, double f_xmin, double f_xmax){
    double x;
    double precision = 0.00000000000001;
    while(fabs((xmax - xmin)/(xmax + xmin)) >= precision){
        x = xmax;
        xmax -= (xmax - xmin) / (f_xmax - f_xmin) * f_xmax;
        xmin = x;
        f_xmin = f_xmax;
        f_xmax = numerov(xmax);
    }
    return (xmax + xmin) / 2;
}

double zeroFinder(double xmin, double xmax, double dx){
    double f_xmin = numerov(xmin);
    double f_xmax;
    for(int i = 1; i <= (int)((xmax-xmin)/dx); i++){
        f_xmax = numerov(xmin + i*dx);
        if(f_xmin * f_xmax < 0){
            return secantMethod(xmin+(i-1)*dx, xmin+i*dx, f_xmin, f_xmax);
        }
        f_xmin = f_xmax;
    }
    printf("Error in zeroFinder: no zero found in the interval (%f, %f) with dx = %f\n", xmin, xmax, dx);
	exit(EXIT_FAILURE);
}

// calculates the value of mu = E + Vint from Numerov solution
double efunc(double * y, double mu){
    double vint1 = pow(y[n],2) * pow(y[n]/(h*n),2) / 2;
    double vint2 = R[n] * pow(y[n]/(h*n),2) / 2;
    for(int i = 1; i < n; i++){
        vint1 += pow(y[i],2) * pow(y[i]/(h*i),2);
        vint2 += R[i] * pow(y[i]/(h*i),2);
    }

    printf("E = %+.15f", mu + Na * h * (vint1/2 - vint2));

    return Na * h * (vint1 - vint2);
}

double step(double * u){
    // update vector to be used in the potential
    for(int i = 0; i <= n; i++){
        R[i] = alpha*pow(u[i],2) + (1-alpha)*R[i];
    }

    // I check the norm with the trapeziudal rule
    /*double norm = (R[0] + R[n]) / 2;
    for(int i = 1; i < n; i++){
        norm += R[i];
    }
    printf("%f\n",sqrt(norm * h));*/

    // compute the new lowest eigenvalue
    double mu = zeroFinder(0, 12, 0.01);

    // update eigenvector
    numerovSave(mu, u);

    // check error on the energy functional
    return efunc(u, mu);
}

int main(){
    //double rmax = h * n; // should be of the order of 10
    double precision = 0.000000000000001;
    R = malloc(sizeof(double)*(n+1)); // to be used in the potential
    double u [n+1]; // solution of the GP equation
    for(int i = 0; i <= n; i++){
        u[i] = 2*(h*i)*exp(-pow(h*i,2)/2)/pow(M_PI,0.25);
        R[i] = pow(u[i], 2);//0;
    }
    double check;
    int counter = 0;
    do{
        counter++;
        check = step(u);
        printf("   (step %03d, error %+.15f)\n", counter, check);
    }while(fabs(check) >= precision);

    return 0;
}




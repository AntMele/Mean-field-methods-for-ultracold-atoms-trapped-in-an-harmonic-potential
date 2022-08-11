
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>


int n = 1000; // i = 0,...,n but since u_0 = u_n = 0 vectors and H will have dimension n-1
double h = 0.01;
double Na = 100;
double alpha = 0.1;
gsl_vector * R; // to be used in the potential


// diagonal H element
double U(int i){
    return 1/pow(h,2) + 0.5*pow(h*i,2) + Na*gsl_vector_get(R,i-1)/pow((h*i),2);
}

double k_square(double E, int i){
    return 2*E - pow(h*i,2) - 2*Na*gsl_vector_get(R,i-1)/pow((h*i),2);
}

// calculates the value of mu = E + Vint from Numerov solution
double efunc(gsl_vector * y, double mu){
    double ekin = 0;
    double vext = 0;
    double vint = 0;
    for(int i = 1; i < n; i++){
        ekin += pow(gsl_vector_get(y,i-1),2) * k_square(mu, i);
        vext += pow(gsl_vector_get(y,i-1)*h*i,2);
        vint += pow(gsl_vector_get(y,i-1),4) / pow(h*i,2);
    }
    ekin /= 2;
    vext /= 2;
    vint *= Na;

   // printf("E = %.15f", mu - h*vint/2);
    printf("mu = %.15f", mu );
    return h * (ekin + vext + vint);
}

double finiteDifferences(gsl_vector * u_new){
    // generate vector to be used in the potential
    for(int i = 0; i < n-1; i++){
        gsl_vector_set(R,i,alpha*pow(gsl_vector_get(u_new,i),2) + (1-alpha)*gsl_vector_get(R,i));
    }


    // declare Hamiltonian matrix and set it to zero
    gsl_matrix * H = gsl_matrix_calloc(n-1,n-1);

    // set first row
    gsl_matrix_set(H,0,0,U(1));
    gsl_matrix_set(H,0,1,-0.5/pow(h,2));

    // set last row
    gsl_matrix_set(H,n-2,n-2,U(n-1));
    gsl_matrix_set(H,n-2,n-3,-0.5/pow(h,2));

    // set the rest of H matrix
    for(int i = 1; i < n-2; i++){
        gsl_matrix_set(H,i,i,U(i+1));
        gsl_matrix_set(H,i,i-1,-0.5/pow(h,2));
        gsl_matrix_set(H,i,i+1,-0.5/pow(h,2));
    }

    // diagonalize H
    gsl_vector * eval = gsl_vector_alloc(n-1);
    gsl_matrix * evec = gsl_matrix_alloc(n-1,n-1);
    gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc(n-1);
    gsl_eigen_symmv(H, eval, evec, w);

    // update eigenvectors
    int j = gsl_vector_min_index(eval);
    gsl_matrix_get_col(u_new, evec, j);
    gsl_vector_scale(u_new, 1/sqrt(h));
    double mu = gsl_vector_get(eval, j);

    return efunc(u_new, mu) - mu;
}

int main(){
    //double rmax = h * n; // should be of the order of 10
    double precision = 0.000000000001;
    R = gsl_vector_alloc(n-1); // to be used in the potential
    gsl_vector * u_new = gsl_vector_calloc(n-1); // solution of the GP equation at the step i (set to zero)
    double check;
    for(int i = 0; i < n-1; i++){
        gsl_vector_set(u_new,i,2*(h*(i+1))*exp(-pow(h*(i+1),2)/2)/pow(M_PI,0.25));
        gsl_vector_set(R,i,pow(gsl_vector_get(u_new, i), 2));
    }

    int counter = 0;
    do{
        counter++;
        check = finiteDifferences(u_new);
        printf("   (step %03d, error %+.15f)\n", counter, check);
    }while(fabs(check) >= precision);

    return 0;
}

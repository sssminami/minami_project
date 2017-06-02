#include <stdio.h>
#include <math.h>

#define step 200

struct object
{
	double x[step];
	double y[step];
	double vx[step];
	double vy[step];
	double ax[step];
	double ay[step];
};


int main(){

	int i,j;
	const double dt = 0.01;
	const double sigma = 3.4e-10;
	const double m = 6.63e-26;
	const double epsilon = 120;
	const double tau_Ar = 0.31e-12;
	const double l = 6.0, n = 12.0;
	double r[step];

	struct object Ar1,Ar2;

	Ar1.x[0] = 0.0;
	Ar1.y[0] = 3.0;
	Ar1.vx[0] = sqrt(3.0);
	Ar1.vy[0] = -3.0;

	Ar2.x[0] = 0.0;
	Ar2.y[0] = -1.0;
	Ar2.vx[0] = sqrt(3.0);
	Ar2.vy[0] = 1.0;

	r[0] = sqrt(pow((Ar1.x[0] - Ar2.x[0]), 2) + pow((Ar1.y[0] - Ar2.y[0]), 2));

	for (i = 1; i <= 2; i++){
		for (j = 1; j < step; j++){
			//calcucation
			Ar1.ax[j] = (pow(1 / r[j - 1], n + 1) - (l / n)*pow(1 / r[j - 1], l + 1))*((Ar1.x[j - 1] - Ar2.x[j - 1]) / r[j - 1]);
			Ar2.ax[j] = (pow(1 / r[j - 1], n + 1) - (l / n)*pow(1 / r[j - 1], l + 1))*((Ar2.x[j - 1] - Ar1.x[j - 1]) / r[j - 1]);
			Ar1.ay[j] = (pow(1 / r[j - 1], n + 1) - (l / n)*pow(1 / r[j - 1], l + 1))*((Ar1.y[j - 1] - Ar2.y[j - 1]) / r[j - 1]);
			Ar2.ay[j] = (pow(1 / r[j - 1], n + 1) - (l / n)*pow(1 / r[j - 1], l + 1))*((Ar2.y[j - 1] - Ar1.y[j - 1]) / r[j - 1]);

			//update
			Ar1.x[j] = Ar1.x[j - 1] + dt*Ar1.vx[j - 1] + (pow(dt, 2) / 2 * m)*m*Ar1.ax[j - 1];
			Ar1.vx[j] = Ar1.vx[j - 1] + (dt / (2 * m))*(m*Ar1.ax[j] + m*Ar1.ax[j - 1]);

			Ar1.y[j] = Ar1.y[j - 1] + dt*Ar1.vy[j - 1] + (pow(dt, 2) / 2 * m)*m*Ar1.ay[j - 1];
			Ar1.vy[j] = Ar1.vy[j - 1] + (dt / (2 * m))*(m*Ar1.ay[j] + m*Ar1.ay[j - 1]);

			Ar2.x[j] = Ar2.x[j - 1] + dt*Ar2.vx[j - 1] + (pow(dt, 2) / 2 * m)*m*Ar2.ax[j - 1];
			Ar2.vx[j] = Ar2.vx[j - 1] + (dt / (2 * m))*(m*Ar2.ax[j] + m*Ar2.ax[j - 1]);

			Ar2.y[j] = Ar2.y[j - 1] + dt*Ar2.vy[j - 1] + (pow(dt, 2) / 2 * m)*m*Ar2.ay[j - 1];
			Ar2.vy[j] = Ar2.vy[j - 1] + (dt / (2 * m))*(m*Ar2.ay[j] + m*Ar2.ay[j - 1]);

			r[j] = sqrt(pow((Ar1.x[j] - Ar2.x[j]), 2) + pow((Ar1.y[j] - Ar2.y[j]), 2));
		}
	}
	printf("---------\n");
	for (j = 0; j < step; j++){
		printf("%lf     %lf     %lf      %lf     %lf\n", j*dt, Ar1.x[j], Ar1.y[j], Ar2.x[j], Ar2.y[j]);
	}
	printf("---------\n");
	for (j = 0; j < step; j++){
		printf("%lf     %lf     %lf      %lf     %lf\n", j*dt, Ar1.vx[j], Ar1.vy[j], Ar2.vx[j], Ar2.vy[j]);
	}


}

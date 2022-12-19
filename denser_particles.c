#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

# define N     20
# define times 30000
# define save  30
# define fps   30 

# define XMIN 0.0
# define XMAX 2.0
# define YMIN 0.0
# define YMAX 2.0

void create_data(char* fname);

void create_gnuplot_file(char* fname, char* data_file);

void create_plots(char* fname, char* data_file,  float margin);

void create_video(char* fname);

int main() {
	
	// Creating the data
	char* data_file = "data.dat";
	create_data(data_file);
	
	// Creating plots
	char* gnuplot_script = "plots.p";
	double margin = 0.25;
	create_plots(gnuplot_script, data_file, margin);
	
	// Createing the video
	char* video_file = "denser_particle_movement";
	create_video(video_file);
	
	printf("Done :D\n");	
}

// ******************************* Functions *******************************
void create_data(char* fname){
		
	printf("*************Writing Data*************\n");
	//    File
	FILE *fp = fopen(fname, "w");

	// ~~~~~~~~~~~~~~~~~~~~~~~ particle much denser than the surrounding fluid ~~~~~~~~~~~~~~~~~~~~~~~
    double pi = 3.141592653589793;
    double St = 0.001;
    double rc = 1;
    double rd = 10;
    double Dp = 0.00006;
	//double Dp = 0.001;
    double rp = Dp / 2;
    double mp = rd * pi / 6 * pow(Dp, 3);
    double mf = 1;
    double Uo = 1;
    double mc = 0.000001;
	//double mc = 0.001;
    double g = 9.81;
    double a = 2;
	//double a = 0.63661977236;
    double R = mf / (mp + mf / 2);
	
    double dx = (XMAX - XMIN) / N;
    double dy = (YMAX - YMIN) / N;

    // Initializing needed matrices
    double x1[N][N] = {0.0};
    double x2[N][N] = {0.0};

    double u1[N][N] = {0.0};
    double u2[N][N] = {0.0};

    double V1[N][N] = {0.0};
    double V2[N][N] = {0.0};

// ~~~~~~~~~~~~~~~~~~~~~~ Equations of Motion ~~~~~~~~~~~~~~~~~~~~~~
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Calculating x1
            x1[i][j] = XMIN + (j + 1) * dx;

            // Calculating X2
            x2[i][j] = YMIN + (i + 1) * dy;

            // Calculating u1 and u2
            u1[i][j] = Uo * cos(x1[i][j] / a) * cos(x2[i][j] / a);
            u2[i][j] = Uo * sin(x1[i][j] / a) * sin(x2[i][j] / a);

            // Calculating V1 and V2
            V1[i][j] = u1[i][j];
            V2[i][j] = u2[i][j];
        }
    }

    double dt = 0.001;

    double Rer1[N][N] = {0.0};
    double Rer2[N][N] = {0.0};

    double f1[N][N] = {0.0};
    double f2[N][N] = {0.0};

    double steady_state_x[N][N] = {0.0};
    double steady_state_y[N][N] = {0.0};
    double external_forces_x[N][N] = {0.0};
    double external_forces_y[N][N] = {0.0};
    double virtual_mass_x[N][N] = {0.0};
    double virtual_mass_y[N][N] = {0.0};

    double dV1[N][N] = {0.0};
    double dV2[N][N] = {0.0};

    double V1_new[N][N] = {0.0};
    double V2_new[N][N] = {0.0};
    double x1_new[N][N] = {0.0};
    double x2_new[N][N] = {0.0};

    double gravity = 0;
    double value1 = 0;
    double value2 = 0;
    double tv = 0;

    double x1_val = -1;
    double x2_val = 0;

    for (int t = 0; t < times; t++) {

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                value1 = abs((u1[i][j] - V1[i][j])) * Dp * rc / mc;
                value2 = abs((u2[i][j] - V2[i][j])) * Dp * rc / mc;

                if (value1 == 0)
                    Rer1[i][j] = pow(10, -7);
                else
                    Rer1[i][j] = value1;
                f1[i][j] = 1;

                if (value2 == 0)
                    Rer2[i][j] = pow(10, -7);
                else
                    Rer2[i][j] = value2;
                f1[i][j] = 1;


                // Calculating dV1
                tv = rd * pow(Dp, 2) / (18 * mc);
                steady_state_x[i][j] = 1 / tv * (Uo * cos(x1[i][j] / a) * cos(x2[i][j] / a) - V1[i][j]);
                external_forces_x[i][j] = -rc / rd * (pow(Uo, 2) / a * cos(x1[i][j] / a) * sin(x1[i][j] / a));
                virtual_mass_x[i][j] = rc / rd * (-V1[i][j] / 2 * Uo / a * sin(x1[i][j] / a) * cos(x2[i][j] / a) -
                                                  V2[i][j] / 2 * Uo / a * cos(x1[i][j] / a) * sin(x2[i][j] / a));
                dV1[i][j] = (1 / (1 + rc / rd)) * (steady_state_x[i][j] + external_forces_x[i][j] + virtual_mass_x[i][j])*dt;


                // Calculating dv2
                gravity = -g;
                steady_state_y[i][j] = 1 / tv * (Uo * sin(x1[i][j] / a) * sin(x2[i][j] / a) - V2[i][j] / 2);
                external_forces_y[i][j] = rc / rd * (pow(Uo, 2) / a * cos(x2[i][j] / a) * sin(x2[i][j] / a) + g);
                virtual_mass_y[i][j] = rc / rd * (V1[i][j] / 2 * Uo / a * cos(x1[i][j] / a) * sin(x2[i][j] / a) +
                                                  V2[i][j] / 2 * Uo / a * sin(x1[i][j] / a) * cos(x2[i][j] / a));

                dV2[i][j] = (1 / (1 + rc / rd)) *
                            (steady_state_y[i][j] + external_forces_y[i][j] + virtual_mass_y[i][j] + gravity) * dt;


                V1_new[i][j] = V1[i][j] + dV1[i][j];
                V2_new[i][j] = V2[i][j] + dV2[i][j];

                x1_new[i][j] = x1[i][j] + (V1[i][j] + V1_new[i][j])*dt/2;
                x2_new[i][j] = x2[i][j] + (V2[i][j] + V2_new[i][j])*dt/2;

                x1_val = x1[i][j] + (V1[i][j] + V1_new[i][j]) * dt / 2;
                x2_val = x2[i][j] + (V2[i][j] + V2_new[i][j]) * dt / 2;

                if (x1_new[i][j] > XMAX)
                    x1_new[i][j] = XMIN + x1_val - XMAX;
                else if (x1_new[i][j] < XMIN)
                    x1_new[i][j] = XMAX - abs(x1_val - XMIN);

                if (x2_new[i][j] > YMAX)
                    x2_new[i][j] = YMIN + x2_new[i][j] - YMAX;
                else if (x2_new[i][j] < YMIN)
                    x2_new[i][j] = YMAX - abs(x2_val - YMIN);
            }
        }

        for (int i=0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                x1[i][j] = x1_new[i][j];
                x2[i][j] = x2_new[i][j];
                V1[i][j] = V1_new[i][j];
                V2[i][j] = V2_new[i][j];
                u1[i][j] = Uo*cos(x1[i][j]/a)*cos(x2[i][j]/a);
                u2[i][j] = Uo*sin(x1[i][j]/a)*sin(x2[i][j]/a);
            }
        }

		if (t % save == 0){
			for (int i=0; i < N; i++) {
				for (int j=0; j < N; j++) {
					fprintf(fp, "%.15f ", x1[i][j]);
					fprintf(fp, "%.15f\n", x2[i][j]);
				}
			}
		}
    }
    fclose(fp);
	printf("Done creating data\n");
}

void create_plots(char* fname, char* data_file, float margin){
	printf("***************Plotting Data***************\n");

	create_gnuplot_file(fname, data_file);
	
	// Creating the plots 	
	char command[80] ={""};
	sprintf(command, "gnuplot -c %s \"%f\" \"%f\" \"%f\" \"%f\" \"%d\" \"%d\" ", fname, XMIN-margin, XMAX+margin, YMIN-margin, YMAX+margin, N*N, times/save);
	system(command);
	
	printf("Done creating the plots\n");
}

void create_gnuplot_file(char* fname, char* data_file){
	FILE* fp = fopen(fname, "w");
	
	fprintf(fp, "set terminal png size 500,500\n");
	fprintf(fp, "set xrange [ARG1:ARG2]\n");
	fprintf(fp, "set yrange [ARG3:ARG4]\n");

	fprintf(fp, "N = ARG5\n");
	fprintf(fp, "times = ARG6\n");

	fprintf(fp, "i = 0\n");
	fprintf(fp, "j = N\n");

	fprintf(fp, "do for [t=1:times]{\n");
	fprintf(fp, "	outfile = sprintf('plots/%c07d.png',t)\n", 37);
	fprintf(fp, "	set output outfile\n");
	fprintf(fp, "	plot \"%s\" every ::i::j notitle\n", data_file);
	fprintf(fp, "	i = i + N\n");
	fprintf(fp, "	j = j + N\n}");	
	
	fclose(fp);
}

void create_video(char* fname){
	printf("***************Videoing Data***************\n");
	
	char command[100] ={""};
	//Creating the video
	sprintf(command, "ffmpeg -loglevel quiet -y -i plots/%c7d.png -filter:v fps=%d %s.mpeg\n\n", 37, fps, fname);
	system(command);
	
	printf("Done creating the video\n");
}
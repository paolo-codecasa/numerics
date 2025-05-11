/*---------------------------------------------*\
|  CWR 2023                                     |
|  Blatt 4 - Aufgabe 14 (PV)                    |
|  Author  p.codecasa@stud.uni-goettingen.de    |
\*---------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "my_numerics.h"
//Es werden alle benötigten Bibliotheken reingeladen.(auch rk4step und eulerstep)

//Definition von ODE_dual springs
void ODE_dual_springs(double t, double *y, double *f, void *params) {
    double xwert = y[0];
    double ywert = y[1];


    double r0 = ((double*)params)[0];
    double k =((double*)params)[1];
    double m =((double*)params)[2];
    double x_A =((double*)params)[3];
    double x_B =((double*)params)[4];
    double y_A =((double*)params)[5];
    double y_B =((double*)params)[6];


    double betrag_A = sqrt((xwert-x_A)*(xwert-x_A)+(ywert-y_A)*(ywert-y_A));
    double betrag_B = sqrt((xwert-x_B)*(xwert-x_B)+(ywert-y_B)*(ywert-y_B));

    double a_x = -(k/m) *( (betrag_A-r0)*((xwert-x_A)/betrag_A) + (betrag_B-r0)*((xwert-x_B)/betrag_B));
    double a_y = -(k/m) *( (betrag_A-r0)*((ywert-y_A)/betrag_A) + (betrag_B-r0)*((ywert-y_B)/betrag_B));


    f[0] = y[2];
    f[1] = y[3];
    f[2] = a_x;
    f[3] = a_y;

    return;
}

//Ab hier geht es los.
int main() {

    double r0 = 1;
    double k = 10;
    double m = 1;
    double x_A = -1;
    double x_B = 1;
    double y_A = 0;
    double y_B = 0;

    double deltat = 1e-4;

    double params[] = {r0,k,m,x_A,x_B,y_A,y_B};
    double y1[] = {-0.5,0,0,0};
    double t0 = 0;
    double t_x = 2*M_PI/sqrt(2*k);
    




//Im folgenden werden alle Daten für die jeweiligen Plots erzeugt

//////////////////////////////////////////////////////////////////////
//4
printf("Begin 4 \n");
//Integration mit euler
    FILE* myFile1 = fopen("schwingereuler.csv", "w");

    while(t0 < 10 * t_x) {

        euler_step(t0,deltat,y1,ODE_dual_springs,4,params);
        fprintf(myFile1,"%g,%g,%g,%g,%g\n",t0,y1[0],y1[1],y1[2],y1[3]);
        t0 += deltat;

    }

    fclose(myFile1);

//////////////////////////////////////////////////////////////////
//5
printf("Begin 5 \n");
//Integration mit rk4
    double y2[] = {-0.5,0,0,0};
    t0 = 0;

    FILE* myFile2 = fopen("schwingerrk4.csv", "w");

    while(t0 < 10 * t_x) {

        rk4_step(t0,deltat,y2,ODE_dual_springs,4,params);
        fprintf(myFile2,"%g,%g,%g,%g,%g\n",t0,y2[0],y2[1],y2[2],y2[3]);
        t0 += deltat;

    }

    fclose(myFile2);

//////////////////////////////////////////////////////////////////
//6
printf("Begin 6 \n");
//Nun wird der Fehler von euler für jeweils logarithmisch verteile kleine delta Schritte errechnet.
    FILE* myFile3 = fopen("schwingereulerlog.csv", "w");

    for(int i = 0; i < 1000;i++) {

        t0 = 0;
        double y3[] = {-0.5,0,0,0};
        deltat = pow(10,-(1 +  (4./999) * i));

        while(t0 < 10 * t_x) {
                //Eine if Bedingung sorgt dafuer dass man immer genau bei 10*t_x den letzten Schritt auswertet.
                if(t0 + deltat > 10*t_x) {
                    deltat = 10*t_x -t0;
                }
            //euler step
            euler_step(t0,deltat,y3,ODE_dual_springs,4,params);
            t0 +=deltat;
        
        }
        //Fehler wird berechnet und mit deltat in Datei geschrieben.
        deltat = pow(10,-(1 +  (4./999) * i));
        double error = fabs(-0.5-y3[0]);
        fprintf(myFile3,"%g,%g\n",deltat,error);
    

    }

    fclose(myFile3);


//Nun wird der Fehler von rk4 für jeweils logarithmisch verteile kleine delta Schritte errechnet.
//////////////////////////////////////////////////////////////////
    FILE* myFile4 = fopen("schwingerrk4log.csv", "w");

    for(int i = 0; i < 1000;i++) {

        t0 = 0;
        double y4[] = {-0.5,0,0,0};
        deltat = pow(10,-(1 +  (4./999) * i));

        while(t0 < 10 * t_x) {
            //Eine if Bedingung sorgt dafuer dass man immer genau bei 10*t_x den letzten Schritt auswertet.
            if(t0 + deltat > 10*t_x) {
                    deltat = 10*t_x -t0;
                }
            //rk4 step
            rk4_step(t0,deltat,y4,ODE_dual_springs,4,params);
            t0 +=deltat;
        
        }
        //Fehler wird berechnet und mit deltat in Datei geschrieben.
        deltat = pow(10,-(1 +  (4./999) * i));
        double error = fabs(-0.5-y4[0]);
        fprintf(myFile4,"%g,%g\n",deltat,error);
    

    }

    fclose(myFile4);

//////////////////////////////////////////////////////////////////
//7
printf("Begin 7 \n");
//Trajektorie1
///////////////////////////////////////////////////////////////////////////


    double y5[] = {-0.5,0.5,0,0};
    t0 = 0;
    FILE* myFile5 = fopen("Trajektorie1.csv", "w");

    while(t0 < 0.5) {

        rk4_step(t0,1e-4,y5,ODE_dual_springs,4,params);
        fprintf(myFile5,"%g,%g\n",y5[0],y5[1]);
        t0 += deltat;

    }

    fclose(myFile5);

//Trajektorie2
//////////////////////////////////////////////////////////////////////////
    double y6[] = {-0.1,0.5,0,0};
    t0 = 0;
    FILE* myFile6 = fopen("Trajektorie2.csv", "w");

    while(t0 < 0.5) {

        rk4_step(t0,1e-4,y6,ODE_dual_springs,4,params);
        fprintf(myFile6,"%g,%g\n",y6[0],y6[1]);
        t0 += deltat;

    }

    fclose(myFile6);

//Trajektorie3
//////////////////////////////////////////////////////////////////////////////
double y7[] = {-0.1,0.1,0,3};
    t0 = 0;
    FILE* myFile7 = fopen("Trajektorie3.csv", "w");

    while(t0 < 0.5) {

        rk4_step(t0,1e-4,y7,ODE_dual_springs,4,params);
        fprintf(myFile7,"%g,%g\n",y7[0],y7[1]);
        t0 += deltat;

    }

    fclose(myFile7);

    return 0;
}



/*
Antworten auf die Fragen aus Aufgabe 14:

Aufgabe 14.5:

    - sowohl der Plott von euler und rk4 ueberlagern sich jeweils ,sodass nur einer Linie zu erkennen ist.
        Es kann also darauf geschlossen werdendass im vorhanden Sachverhaltbei den verwendeten Parametern beide integratoren sehr gut sind.
      bei der Inegration mit euler ist bei spaeteren Zeiten, ab sekunde 8 bei den Minima der Schwingun ein kleiner unterschied zwischen
      grüner(numerisch) und analytischer(rot) kurve zu sehen. Beim rk4 verfahren jedoch nicht. ==> In dieser konstellation rk4 kaum merklich aber etwas besser.

Aufgabe14.6

    - Es ist zu erkennen das beide Integratoren für kleinere deltas viel genauer am analytischen ergebnis liegen.
      Es fällt jedoch auf, dass rk4 sich schon früher (schon bei größeren deltas) dem analytischen Ergebnis naehert und für sehr kleine deltas die rk4 fehler immer kleiner als die eulerfehler sind.(Rote gerade unter grüner gerade).
      Die Fehlergerade von rk4 fällt außerdem schneller ab als die Fehlergerade vin euler.(Die grüne Gerade hat eine höhere Steigung als die rote)
      Außerdem sind im deltabeecih von 1e-1 bis 1e-2 bei der eulerintegration verauschte punkte zu erkennen, hier sind die deltas also noch zu groß sodass euler sehr variert. Bei rk4 ist dies jedoch schon nicht mehr der Fall.
      Des Weiteren ist beim rk4 im delta Bereich von 1e-3 bis 1e-5 auch ein rauschen zu erkennen, dies könnte an Unterlauefen liegen, da die Fehler sehr sehr gering (1e-15 Bereich) sind.
      ==> rk4 viel besser als euler.

Aufgabe14.7

    - Nach Aufgabe14.6 ist der Integrator rk4 besser und wird daher verwendet fuer die Simulation von Trajektorien (bis halbe sekunde simuliert).

      Trajektorie1(x = -0.5, y = 0.5, v_x = 0, v_y = 0):

        - Die Massekugel schwingt in der x Richtung zwischen -0.5 und 0.5 hin und her. Die y _Komponente bleibt jedoch immer positiv.
          Es ist eine schöner geschwungende Figur zu erkennen die wie ein nach links offenes U ausseiht.

      Trajektorie2(x = -0.1, y = 0.5, v_x = 0, v_y = 0):

        - Die Massekugel schwingt in der x Richtung zwischen -0.1 und 0.1 hin und her. Die y _Komponente geht dabei von o.5 bis fast -0.5 und dann wieder zurueck zu 0.5.
          Es ist eine Figur zu erkennen in der Elipsen eingeschlossen sind. Kreuzen der Lininen inder Naehe des Ortes x = 0.

    Trajektorie3(x = -0.1, y = 0.1, v_x = 0, v_y = 3):

        - Die Massekugel schwingt vom Anfangsort zunaechst nach oben und in positive x- Richtung, da der Massepunkt eine Startgeschwindigkeit in positive y -Richtung besitz.
          Anschliessend kehrt der Massepunkt um und geht in negative y und negative x Richtung. Der wweitere Verlauf ist scher zu beschreiben, da nur noch schwer ein Muster erkennbar ist. Schon etwas wirrer als die vorherigen Plots.

*/


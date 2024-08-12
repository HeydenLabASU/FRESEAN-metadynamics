#include <stdio.h>

int progressBar(int current,int interval,int total) {
        int p1,perCent;
        int i;

        p1=current+1;
        if(p1%interval==0) {
                perCent=(100*p1)/total;
                printf("\r");
                i=0;
                while(i<perCent) {
                        printf("*");
                        i+=2;
                }
		while(i<100) {
                        printf(" ");
                        i+=2;
                }
                printf("[%3d%%]",perCent);fflush(stdout);
        }
        if(p1==total) {
                printf("\n");
                fflush(stdout);
        }
}

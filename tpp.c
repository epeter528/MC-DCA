#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>




int main(int argc, char *argv[]) 

{
    
    int i,j,k,l,m,n;
    FILE *fp2;
    char string[100];
    char fpname[100];
    int  *x1,*y1;
    int  *x2,*y2;
    int  x,y;
    int size = 1000000;
    int number_of_contacts;
    int number_of_real_contacts;
    float real_contact;
    
    x1 = (int *) malloc ( sizeof(int)*size);
    y1 = (int *) malloc ( sizeof(int)*size);
    x2 = (int *) malloc ( sizeof(int)*size);
    y2 = (int *) malloc ( sizeof(int)*size);    
    
  //  printf("%s\t%s\n",argv[0],argv[1]);
    
    strcpy(fpname,argv[1]);
    
    FILE *fp = fopen(fpname,"r");

    i = 1;
    
    while(fgets(string,100,fp)!=0) {
        
        sscanf(string,"%d%d",&x,&y);
        
        x1[i] = x;
        y1[i] = y;
        
    //    printf("%d\t%d\t%d\n",x1[i],y1[i],i);
        
        i++;
    };
    
    fclose(fp);
    
    number_of_real_contacts = i-1;
    
    strcpy(fpname,argv[2]);
    
    fp = fopen(fpname,"r");

    i = 1;
    
    while(fgets(string,100,fp)!=0) {
        
        sscanf(string,"%d%d",&x,&y);
        
        x2[i] = x;
        y2[i] = y;
        
      //  printf("%d\t%d\t%d\n",x2[i],y2[i],i);
        
        i++;
    };
    
    fclose(fp);

    number_of_contacts = i-1;

  //  printf("%d\t%d\n",number_of_contacts,number_of_real_contacts);  

    real_contact = 0.0;
    
    for(k=1;k<=number_of_real_contacts;k++) {
    
        for(i=1;i<=number_of_contacts;i++) {
        
            if( x1[k] - x2[i] == 0 && y1[k] - y2[i] == 0) {
                real_contact = real_contact + 1.0;
                
                break;
            
            };
 //           if( x1[i] - y2[k] == 0 && y1[i] - x2[k] == 0) {
 //               real_contact = real_contact + 1.0;
            
 //               break;
 //           };            
            
        };
        
        printf("%d\t%f\n",k,real_contact/((float)k));
        
    };
};

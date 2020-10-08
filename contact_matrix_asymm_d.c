#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>




int main(int argc, char *argv[]) 

{
    
    int i,j,k,l,m,n;
    char string[100];
    char fpname[100];
    char stringA[100];
    int resnum_A,resnum_B;
    char restype_A[2],restype_B[2];
    float MI_score,corrected_norm_score;
    int *resnum_A_arr,*resnum_B_arr;
    float *MI_score_arr,*corr_MI_score_arr;
    int resnumA_max=0;
    int number_of_elements=0;
    int threshold=0;
    int size;
    int counter = 0;
    int dist=4;

    FILE *fp2;
    
    size = 10000;
    i    = 1;
    
    fp2 = fopen("contact_corr_MI_2d.xvg","w");
    
    printf("%s\n","number of considered contacts ");
    scanf("%d", &threshold);
    printf("%s\n","Distance threshold between residues ");
    scanf("%d", &dist);


    resnum_A_arr = malloc((sizeof(int))*size);
    resnum_B_arr = malloc((sizeof(int))*size);
    MI_score_arr = malloc((sizeof(float))*size);
    corr_MI_score_arr = malloc((sizeof(float))*size);
    
    printf("%s\t%s\n",argv[0],argv[1]);
    
    strcpy(fpname,argv[1]);

    FILE *fp = fopen(fpname,"r");

    n = 1;
    
    while(fgets(string,100,fp)!=0) {
    
    strcpy(stringA,string);
            
    printf("%s\t%d\n",stringA,i);
            
   if(i > 0) {         
        
        sscanf(string,"%d%d%f",&resnum_A,&resnum_B,&MI_score);
    //    printf("%d\t%d\t%f\n",resnum_A,resnum_B,MI_score);
    
     if(n <= threshold) {   
        
      if(resnum_A - resnum_B >= dist || resnum_A - resnum_B <= -dist) {  
          
        resnum_A_arr[n] = resnum_A;
        resnum_B_arr[n] = resnum_B;
        MI_score_arr[n] = MI_score;  
        
        fprintf(fp2,"%d\t%d\n",resnum_A,resnum_B);        
        n++;
        
      };
     };  
      
        printf("%d\n",n);
        
    //   return 0;
        
        if(resnum_A > resnumA_max) {
       
           resnumA_max = resnum_A;

          };
         };
   
        i++;
        
    };    
    
    fclose(fp);
    fclose(fp2);
       
}

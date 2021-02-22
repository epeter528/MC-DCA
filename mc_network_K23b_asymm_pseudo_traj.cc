#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <locale>
#include <omp.h>
#include <random>

int partition(double *a,int *b, int *c,int start,int end)
{
    double pivot=a[end];
    //P-index indicates the pivot value index
    
    int P_index=start;
    int i; //t is temporary variable
    
    double t;
    
    int buff_a, buff_b;
    
    //Here we will check if array value is 
    //less than pivot
    //then we will place it at left side
    //by swapping 
    
for(i=start;i<end;i++)
{
    	if(a[i]>=pivot)
        {
            buff_a = b[i];
            b[i] = b[P_index];
            b[P_index] = buff_a;
            
            buff_b = c[i];
            c[i] = c[P_index];
            c[P_index] = buff_b;
            
            t=a[i];
            a[i]=a[P_index];
            a[P_index]=t;
            P_index++;
        }
     }
     //Now exchanging value of
     //pivot and P-index
      buff_a = b[end];
      b[end] = b[P_index];
      b[P_index] = buff_a;
      
      buff_b = c[end];
      c[end] = c[P_index];
      c[P_index] = buff_b;      
     
      t=a[end];
      a[end]=a[P_index];
      a[P_index]=t;
    
     //at last returning the pivot value index
     return P_index;
}
void Quicksort(double *a,int *b, int *c,int start,int end)
{
    if(start<end)
    {
         int P_index=partition(a,b,c,start,end);
             Quicksort(a,b,c,start,P_index-1);
             Quicksort(a,b,c,P_index+1,end);
    }
}

int partition2(double *a,int *b, int *c, int *d,int start,int end)
{
    double pivot=a[end];
    //P-index indicates the pivot value index
    
    int P_index=start;
    int i; //t is temporary variable
    
    double t;
    
    int buff_a, buff_b;
    
    //Here we will check if array value is 
    //less than pivot
    //then we will place it at left side
    //by swapping 
    
for(i=start;i<end;i++)
{
    	if(a[i]>=pivot)
        {
            buff_a = b[i];
            b[i] = b[P_index];
            b[P_index] = buff_a;
            
            buff_b = c[i];
            c[i] = c[P_index];
            c[P_index] = buff_b;
            
            buff_a = d[i];
            b[i] = d[P_index];
            d[P_index] = buff_a;            
            
            t=a[i];
            a[i]=a[P_index];
            a[P_index]=t;
            P_index++;
        }
     }
     //Now exchanging value of
     //pivot and P-index
      buff_a = b[end];
      b[end] = b[P_index];
      b[P_index] = buff_a;
      
      buff_b = c[end];
      c[end] = c[P_index];
      c[P_index] = buff_b;  
      
      buff_b = d[end];
      d[end] = d[P_index];
      d[P_index] = buff_b;
     
      t=a[end];
      a[end]=a[P_index];
      a[P_index]=t;
    
     //at last returning the pivot value index
     return P_index;
}
void Quicksort2(double *a,int *b, int *c, int *d,int start,int end)
{
    if(start<end)
    {
         int P_index=partition2(a,b,c,d,start,end);
             Quicksort2(a,b,c,d,start,P_index-1);
             Quicksort2(a,b,c,d,P_index+1,end);
    }
}


int main(int argc, char **argv) {

    
    int n,i,j,k,l;
    int m,n3;
    
    if (argc < 2) {
        std::cerr << " Wrong format: " << argv[0] << " [infile] " << std::endl;
        return -1;
    }

    std::ifstream input(argv[1]);
    if (!input.good()) {
        std::cerr << "Error opening: " << argv[1] << " . You have failed." << std::endl;
        return -1;
    }
    std::string line, id, DNA_sequence;
    std::string String;
    std::vector <std::string > Sequence_vector;
    char  inputstring[100];
    double tolerance=1E-6;
    double **msa_vector;
    int length;
    double ***msa_vec_d;
    double **single_frequency_count;
    double ***double_frequency_count;
  //  double double_frequency_count[150][150][6][6]; 
    double *weight;
    double seq_id=0.8;
    double q_parameter;
    double pseudocount = 0.5;
    double *G_of_r,**h,**state,****J,****state_2;
    double **f1,****f12;
    double **delta_h,****delta_J;
    double **logPot;
    int    n2,n1;
    int    maxcontact = 1000;
    q_parameter = 5.0;
    double norm;
    int number_of_min_steps = 200;
    int num_iterations = 2;
    double **d_G_dh,****d_G_dJ;
    int    switchflag=1;
    double stepsize_H = 1E-2;
    double stepsize_J = 1E-2;
    double threshold_B=0.95;
    double temperature=1.0;  
    int q_int = 5;
    int number_mc = 10000;
    double lambda=0.01;
    int reseed_interval;
    
    double reg_J, reg_H;

    FILE *fp, *fp2;

    n = 0;
    
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);    
    
#if defined(_OPENMP)    
#pragma omp parallel 
#endif       
    
    
#pragma omp critical    
    if(argc > 1) {
        
    strcpy(inputstring,argv[2]);    
    if(strcmp(inputstring,"rna") == 0) switchflag = 1;
    if(strcmp(inputstring,"protein") == 0) switchflag = 2;
    strcpy(inputstring,argv[3]);
    sscanf(inputstring,"%lg",&temperature);
    strcpy(inputstring,argv[4]);
    sscanf(inputstring,"%d",&number_mc);
    strcpy(inputstring,argv[5]);
    sscanf(inputstring,"%d",&reseed_interval);    
    strcpy(inputstring,argv[6]);
    sscanf(inputstring,"%lg",&stepsize_H);
    strcpy(inputstring,argv[7]);
    sscanf(inputstring,"%lg",&stepsize_J); 
    
   };

   if(switchflag == 2) {
       q_int = 21;
       q_parameter = 21.0;  
   };   
    
    // Parse MSA, argv[1]
    
      while (std::getline(input, line)) {

            id = line.substr(1);
            DNA_sequence.clear();
            DNA_sequence += line;     
            
         if(line[0] != '>') { 
            if(!id.empty())
                std::cout << id << " : " << DNA_sequence << std::endl;
            if(!id.empty()) {
                std::cout << DNA_sequence.length() << " seq. length " << std::endl;
                Sequence_vector.push_back(DNA_sequence);
                length = DNA_sequence.length();
                n++;
            };
         }; 
      };
   
    fp = fopen("DI_ij.dat","w");
    fp2 = fopen("DI_ijk.dat","w");

    // Initialize MSA Vector
    
      std::cout << n << " Number of lines " << std::endl;
      
      int number_of_sequences = n-1;
      
      n = 0;
      while(n <= number_of_sequences) {
        
          std::cout << Sequence_vector[n] << std::endl;
          n++;
          
      };
      
      std::cout << length << " length  " << number_of_sequences << " number of sequences "<< std::endl;
      
      double *excluded;
      
      excluded = (double *) malloc(sizeof(double)*(number_of_sequences+1));
      
      msa_vector = (double **) malloc(sizeof(double)*(number_of_sequences*2));
      
      for(n=0;n<=number_of_sequences+1;n++) {
        
          msa_vector[n] = (double *) malloc(sizeof(double)*(length+1));
          
      };
      
      msa_vec_d = (double ***) malloc(sizeof(double)*(number_of_sequences*2));
      
      for(n=0;n<=number_of_sequences+1;n++) {
        
          msa_vec_d[n] = (double **) malloc(sizeof(double)*(length+1));
          
          for(n2=0;n2<=length;n2++) {
            
              msa_vec_d[n][n2] = (double *) malloc(sizeof(double)*(q_int+1));
              
          };
      };      
      
      for(n=0;n<=number_of_sequences*2-1;n++) {
          for(i=0;i<String.length();i++) {
            
              msa_vector[n][i] = 0;
              
              for(n2=0;n2<=q_int;n2++) {
                
                  msa_vec_d[n][i][n2] = 0.0;
                  
              };
          };
      };
      
     // Set nucleotide indices for each entry in each sequence in a matrix Nxi
      
      n = 0;
      
      while(n <= number_of_sequences) {
          
          String.clear();
          String = Sequence_vector[n];
          
          for(i=0;i<String.length();i++) {
              
         //     std::cout << String.length() << " i " << i << " n " << n << " str.len " << std::endl;
              
            if(switchflag == 1) {  
              
              if(String[i] == 'A') {
                  msa_vector[n][i] = 1;
                  msa_vec_d[n][i][1] = 1.0;
              };
              if(String[i] == 'G') {
                  msa_vector[n][i] = 2;
                  msa_vec_d[n][i][2] = 1.0;                  
              };
              if(String[i] == 'C') {
                  msa_vector[n][i] = 3;
                  msa_vec_d[n][i][3] = 1.0;                  
              };
              if(String[i] == 'U') {
                  msa_vector[n][i] = 4;
                  msa_vec_d[n][i][4] = 1.0;                  
              };
              if(String[i] == '-') {
                  msa_vector[n][i] = 5; 
                  msa_vec_d[n][i][5] = 1.0;                  
              };
         //     std::cout << msa_vector[n+1][i];
            };
            
             if(switchflag == 2) {  
              
              if(String[i] == 'A') {
                  
                  msa_vector[n][i] = 1;
                  msa_vec_d[n][i][1] = 1.0;
              };
              if(String[i] == 'C') {
                  msa_vector[n][i] = 2;
                  msa_vec_d[n][i][2] = 1.0;
              }
              if(String[i] == 'D') {
                  msa_vector[n][i] = 3;
                  msa_vec_d[n][i][3] = 1.0;
              }    
              if(String[i] == 'E') {
                  msa_vector[n][i] = 4;
                  msa_vec_d[n][i][4] = 1.0;
              }    
              if(String[i] == 'F') {
                  msa_vector[n][i] = 5; 
                  msa_vec_d[n][i][5] = 1.0;
              }    
              if(String[i] == 'G') {
                  msa_vector[n][i] = 6;
                  msa_vec_d[n][i][6] = 1.0;
              }   
              if(String[i] == 'H') {
                  msa_vector[n][i] = 7;
                  msa_vec_d[n][i][7] = 1.0;
              }    
              if(String[i] == 'I') {
                  msa_vector[n][i] = 8;
                  msa_vec_d[n][i][8] = 1.0;
              }   
              if(String[i] == 'K') {
                  msa_vector[n][i] = 9;
                  msa_vec_d[n][i][9] = 1.0;
              }    
              if(String[i] == 'L') {
                  msa_vector[n][i] = 10; 
                  msa_vec_d[n][i][10] = 1.0;
              }    
              if(String[i] == 'M') {
                  msa_vector[n][i] = 11;
                  msa_vec_d[n][i][11] = 1.0;
              }    
              if(String[i] == 'N') {
                  msa_vector[n][i] = 12;
                  msa_vec_d[n][i][12] = 1.0;
              }    
              if(String[i] == 'P') {
                  msa_vector[n][i] = 13;
                  msa_vec_d[n][i][13] = 1.0;
              }    
              if(String[i] == 'Q') {
                  msa_vector[n][i] = 14;
                  msa_vec_d[n][i][14] = 1.0;
              }   
              if(String[i] == 'R') {
                  msa_vector[n][i] = 15; 
                  msa_vec_d[n][i][15] = 1.0;                  
              }    
              if(String[i] == 'S') {
                  msa_vector[n][i] = 16;
                  msa_vec_d[n][i][16] = 1.0;                  
              }    
              if(String[i] == 'T') {
                  msa_vector[n][i] = 17;
                  msa_vec_d[n][i][17] = 1.0;
              }    
              if(String[i] == 'V') {
                  msa_vector[n][i] = 18;
                  msa_vec_d[n][i][18] = 1.0;
              }    
              if(String[i] == 'W') {
                  msa_vector[n][i] = 19;
                  msa_vec_d[n][i][19] = 1.0;
              }    
              if(String[i] == 'Y') {
                  msa_vector[n][i] = 20; 
                  msa_vec_d[n][i][20] = 1.0;
              }    
              if(String[i] == '-') {
                  msa_vector[n][i] = 21; 
                  msa_vec_d[n][i][21] = 1.0;
               }   
         //     std::cout << msa_vector[n+1][i];
            };           
              
          };
          
        //  std::cout << "\n";
          
          n++;
      };
      
      weight = (double *) malloc(sizeof(double)*number_of_sequences+1);
      
      // Calculate sequence weights 
      
      double id2,id3;
      
      std::cout << " Calculate sequence weights over " << number_of_sequences << " sequences. Threshold :: "<< seq_id << std::endl;
      
      for(i=0;i<number_of_sequences;i++) {
        
          weight[i] = 0.0;
          excluded[i] = 0.0;
          
      };
    
      
      for(i=0;i<number_of_sequences;i++) {
          
          id2 = 0.0;
          weight[i] = 0.0;
        //  excluded[i] = 0.0;
    
        //  #pragma omp parallel for schedule(dynamic)
          
        for(j=0;j<number_of_sequences;j++) {
          
          if( i != j) {    
              
           for(k=0;k<length;k++) {
               
               if(msa_vector[i][k] == msa_vector[j][k]) {
                   
                   id2 = id2 + 1.0;       
                   
                  };
                 };
             };
           };
       
        weight[i] = id2/((double)(length*number_of_sequences));
       
       };
    //  Sequence_vector.clear();
     

double B_eff,mean,max = 0.0,min=1000.0;
double *m_a;

m_a = (double *) malloc(sizeof(double)*number_of_sequences);

B_eff = 0.0;
mean  = 0.0;
     
for(i=0;i<number_of_sequences;i++) {

 if(weight[i] <= 0.8) m_a[i] = 1.0;	

  for(j=0;j<number_of_sequences;j++) {	
 
    if(i != j) {

      if(weight[j] > 0.8 && weight[i] > 0.8) {
      
	 m_a[i] += 1.0;    

      }
    }
  }
}

for(i=0;i<number_of_sequences;i++) {

    weight[i] = 1.0/m_a[i];

    B_eff     += weight[i];
          
};
     
double **delta_E_N; 

std::cout << "Average over sequences " << std::endl;

J = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  J[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    J[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          J[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};      

double **av_h,****av_J,**corr_h,**d_corr_h,****d_corr_J,****corr_J,**stamp_h,****stamp_J;

h = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    h[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

av_h = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    av_h[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

corr_h = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    corr_h[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

d_corr_h = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    d_corr_h[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

stamp_h = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    stamp_h[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

stamp_J = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  stamp_J[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    stamp_J[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          stamp_J[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
}; 

corr_J = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  corr_J[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    corr_J[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          corr_J[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
}; 

d_corr_J = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  d_corr_J[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    d_corr_J[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          d_corr_J[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
}; 

av_J = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  av_J[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    av_J[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          av_J[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};   

delta_E_N = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    delta_E_N[l] = (double *) malloc(sizeof(double)*(q_int+1));

};


f12 = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  f12[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    f12[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          f12[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};      

f1 = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    f1[l] = (double *) malloc(sizeof(double)*(q_int+1));

};


delta_J = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  delta_J[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    delta_J[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          delta_J[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};      

delta_h = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    delta_h[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

int select1a[length][q_int+1];
int select1b[length][q_int+1];
int select2b[length][q_int+1];


double ****delta_J2, **delta_h2;
double **shift1,****shift12,**prob;

shift12 = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  shift12[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    shift12[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          shift12[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};      

shift1 = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    shift1[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

prob = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    prob[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

delta_J2 = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  delta_J2[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    delta_J2[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          delta_J2[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};    

double ****delta_J3;

delta_J3 = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  delta_J3[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    delta_J3[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          delta_J3[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};   

delta_h2 = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    delta_h2[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

state_2 = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  state_2[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
    state_2[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          state_2[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};   

double **accept_1;

accept_1 = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    accept_1[l] = (double *) malloc(sizeof(double)*(q_int+1));

};

double ****accept_2;

accept_2 = (double ****) malloc(sizeof(double)*(length+1));
    
for(k=0;k<length;k++) {
        
  accept_2[k] = (double ***) malloc(sizeof(double)*(length+1));
         
  for(l=0;l<length;l++) {
           
     accept_2[k][l] = (double **) malloc(sizeof(double)*(q_int+1));
             
     for(n=1;n<=q_int;n++) {             
        
          accept_2[k][l][n] = (double *) malloc(sizeof(double)*(q_int+1)); 
                
        };   
    };
};    


state = (double **) malloc(sizeof(double)*(length+1));

for(l=0;l<length;l++) {
           
    state[l] = (double *) malloc(sizeof(double)*(q_int+1));

};


for(k=0;k<length;k++) {  
    
    for(n=1;n<=q_int;n++) { 

      f1[k][n] = 0.5/((double)q_int)*1.0/(B_eff + 0.5);  
      accept_1[k][n] = 1.0;
      av_h[k][n] = 0.0;
      corr_h[k][n] = 0.0;
      d_corr_h[k][n] = 0.0;      
        
      for(l=0;l<length;l++) {  
        
        for(n2=1;n2<=q_int;n2++) {        
   
            f12[k][l][n][n2] = 0.5/pow((double)q_int,2)*1.0/(B_eff + 0.5);
            accept_2[k][l][n][n2] = 1.0;
            av_J[k][l][n][n2] = 0.0;
            corr_J[k][l][n][n2] = 0.0;
            d_corr_J[k][l][n][n2] = 0.0;
            
        };
      };
    };
};
            
for(i=1;i<number_of_sequences;i++) {  
    
     for(k=0;k<length;k++) {  
    
       for(n=1;n<=q_int;n++) { 

         if(msa_vec_d[i][k][n] != msa_vec_d[i-1][k][n]) f1[k][n] += weight[i]/(B_eff + 0.5);  
          
         for(l=0;l<length;l++) {  
        
           for(n2=1;n2<=q_int;n2++) {
                    
                if(msa_vec_d[i][k][n] != msa_vec_d[i-1][k][n] && msa_vec_d[i][l][n2] != msa_vec_d[i-1][l][n2]) f12[k][l][n][n2] += weight[i]/pow(B_eff + 0.5,1);  
                

             };
           };
         };   
       };  
};

double ran;

for(k=0;k<length;k++) { 
    
   for(n=1;n<=q_int;n++) {
       
      h[k][n]    = f1[k][n];
      stamp_h[k][n] = f1[k][n];

      for(l=0;l<length;l++) {  
        
         for(n2=1;n2<=q_int;n2++) {    
         
           J[k][l][n][n2] = f12[k][l][n][n2];
           stamp_J[k][l][n][n2] = f12[k][l][n][n2];  
           
        }; 
       };
      };
     };

double par_corr_1;
double par_f1;
double par_corr_12;
double par_f12;

double renorm,buff,buff2;

par_f1 = 1.0;
par_f12 = 1.0;
par_corr_1 = 1.0E+0;
par_corr_12 = 1.0E+0;

double Z1,Z2;
double count_reseed = 1.0;

//Z1 = (double *) malloc(sizeof(double)*(length+1));
//Z2 = (double *) malloc(sizeof(double)*(length+1));

double ran2,ran3,*delta_E,**energy1,**energy2;
 
//stepsize_H /= (double)(number_mc);
//stepsize_J /= (double)(number_mc);

double **a1,**a2;

a1 = (double **) malloc(sizeof(double)*(length+1));

for(k=0;k<length;k++) {
  
     a1[k] = (double *) malloc(sizeof(double)*(q_int+1));

}

a2 = (double **) malloc(sizeof(double)*(length+1));

for(k=0;k<length;k++) {

     a2[k] = (double *) malloc(sizeof(double)*(q_int+1));

}

    double *av_k2,*av_l2;
    
    double *av_k3,*av_l3,*av_m3;
    
    av_k2 = (double *) malloc(sizeof(double)*(length+1));
    av_l2 = (double *) malloc(sizeof(double)*(length+1));


    av_k3 = (double *) malloc(sizeof(double)*(length+1));
    av_l3 = (double *) malloc(sizeof(double)*(length+1));    
    av_m3 = (double *) malloc(sizeof(double)*(length+1));  
    
double ***V_coupling;
     
V_coupling = (double ***) malloc(sizeof(double)*(length+1));
     
for(k=0;k<length;k++) {
        
         V_coupling[k] = (double **) malloc(sizeof(double)*(length+1));
         
         for(m=0;m<length;m++) {     
       
            V_coupling[k][m] = (double *) malloc(sizeof(double)*(length+1)); 
             
         };
        };

double **J_coupling;
     
J_coupling = (double **) malloc(sizeof(double)*(length+1));
     
for(k=0;k<length;k++) {
        
         J_coupling[k] = (double *) malloc(sizeof(double)*(length+1));
         
};
     
double ***av_k,***av_l;
double **av_kl;
    
    av_k = (double ***) malloc(sizeof(double)*(length+1));
    
     for(k=0;k<length;k++) {
        
         av_k[k] = (double **) malloc(sizeof(double)*(length+1));
         
         for(l=0;l<length;l++) {
           
             av_k[k][l] = (double *) malloc(sizeof(double)*(q_int+1));
             
         };
      };     
    
    av_l = (double ***) malloc(sizeof(double)*(length+1));
    
     for(k=0;k<length;k++) {
        
         av_l[k] = (double **) malloc(sizeof(double)*(length+1));
         
         for(l=0;l<length;l++) {
           
             av_l[k][l] = (double *) malloc(sizeof(double)*(q_int+1));
             
         };
      };    
    
    av_kl = (double **) malloc(sizeof(double)*(length+1));
     
     for(k=0;k<length;k++) {
        
         av_kl[k] = (double *) malloc(sizeof(double)*(length+1));
         
      };  
      
    double *av_m2;
    
    av_m2 = (double *) malloc(sizeof(double)*(length+1));  
    

for (n1 = 1; n1 <= number_mc ; n1 ++) {   

if(n1 == 1) {      	
	
  for(k=0;k<length;k++) {  
   for(n=1;n<=q_int;n++) {
       
    accept_1[k][n] = 0.0;
       
    for(l=0;l<length;l++) {  
     for(n2=1;n2<=q_int;n2++) {  
         
         accept_2[k][l][n][n2] = 0.0;
         
     };
    };
   };
  };  
};    
    
  Z1 = 0.0;
  Z2 = 0.0;  
  
  for(k=0;k<length;k++) {  
      
   for(n=1;n<=q_int;n++) {  
       
      ran = distribution(generator);    

   //   state[k][n] = ((ran - 0.5)*(h[k][n] - f1[k][n] + 1.0/(0.5+B_eff)))*stepsize_H;

      state[k][n] = f1[k][n] + ((ran-0.5)*f1[k][n])*stepsize_H;

      prob[k][n]      = h[k][n]; 
      
      delta_h2[k][n]  = state[k][n];
      
      for(l=0;l<length;l++) {  
        
         for(n2=1;n2<=q_int;n2++) {      
         
          ran = distribution(generator);

    //      state_2[k][l][n][n2] = ((ran - 0.5)*(J[k][l][n][n2] - f12[k][l][n][n2] + 1.0/(0.5+B_eff)) - 
    //		  (ran - 0.5)*((h[k][n] - f1[k][n] + 1.0/(0.5+B_eff))*(h[l][n2] - f1[l][n2] + 1.0/(0.5+B_eff))))*stepsize_J;

          state_2[k][l][n][n2] = f12[k][l][n][n2] + (ran-0.5)*(f12[k][l][n][n2])*stepsize_J;

          delta_J2[k][l][n][n2] = (J[k][l][n][n2] + state_2[k][l][n][n2]);
          
          prob[k][n]           += J[k][l][n][n2];
          
          delta_h2[k][n]       += state_2[k][l][n][n2];
          
        }; 
       };
       
       Z1 += exp((prob[k][n]));
       Z2 += exp((delta_h2[k][n]));
       
      };
    };
    
//#pragma omp parallel for schedule(dynamic)     
   
  double delta_prob,Z3,Z4;  
  
  Z3 = 0.0;
  Z4 = 0.0;

  for(k=0;k<length;k++) {

   for(n=1;n<=q_int;n++) {

      a1[k][n] = (sqrt(pow(1.0/Z1*(exp(prob[k][n])) - f1[k][n],2)));
      a2[k][n] = (sqrt(pow(1.0/Z2*(exp(delta_h2[k][n])) - f1[k][n],2)));

       for(l=0;l<length;l++) {

         for(n2=1;n2<=q_int;n2++) {

             a1[k][n] += sqrt(pow(1.0/(Z1)*(exp(prob[l][n2])*exp(prob[k][n])) + 1.0/pow(Z1,2)*exp(prob[l][n2])*exp(prob[k][n]) - f12[k][l][n][n2],2));
             a2[k][n] += sqrt(pow(1.0/(Z2)*(exp(delta_h2[l][n2])*exp(delta_h2[k][n])) + 1.0/pow(Z2,2)*exp(delta_h2[l][n2])*exp(delta_h2[k][n]) - f12[k][l][n][n2],2));
	     }
       }

    Z3 += exp(-a1[k][n]/temperature);
    Z4 += exp(-a2[k][n]/temperature);

   }
  }
  
  for(k=0;k<length;k++) {  

   for(n=1;n<=q_int;n++) {

   // printf("%lg\n",1.0/(Z2*temperature)*exp(state[k][n]/temperature));
    
   // Z2 = Z1 - exp(prob[k][n]) + exp(delta_h2[k][n]); 
      
    delta_prob = exp(-(sqrt(pow(1.0/(Z2)*
                       exp(delta_h2[k][n]) - f1[k][n],2))
                     - sqrt(pow(1.0/Z1*exp(prob[k][n]) - f1[k][n],2)))/temperature)*Z3/Z4*exp(((state[k][n] - h[k][n])/temperature));
       
    ran = distribution(generator); 
    
    if(ran <= delta_prob) { //- reg_H*h[k][n]/accept_1[k][n])) {
             
             h[k][n] = state[k][n];
         
             accept_1[k][n] += 1.0;
         
       for(l=0;l<length;l++) {  
        
         for(n2=1;n2<=q_int;n2++) { 
           
           delta_prob = exp(-(sqrt(pow(1.0/(Z2)*(exp(delta_h2[l][n2])*exp(delta_h2[k][n])) + 1.0/pow(Z2,2)*exp(delta_h2[l][n2])*exp(delta_h2[k][n]) - f12[k][l][n][n2],2)) - 
           sqrt(pow(1.0/(Z1)*(exp(prob[l][n2])*exp(prob[k][n])) + 1.0/pow(Z1,2)*exp(prob[l][n2])*exp(prob[k][n]) - f12[k][l][n][n2],2)))/temperature)*Z3/Z4*exp(((state_2[k][l][n][n2] - J[k][l][n][n2]))/temperature);

           ran = distribution(generator); 
           
           if(ran <= delta_prob)  { 
                
                J[k][l][n][n2] = state_2[k][l][n][n2];
              //  h[k][n]        += state[k][n]*state[l][n2];
	      //  h[l][n2]       += state[k][n]*state[l][n2];

                accept_2[k][l][n][n2] += 1.0;
          
         };
        };
       };
     };
    };
   };
 
if(n1%1000 == 0) 

 {

     std::cout << n1 << " mc-steps " << std::endl; 

double d_a1,d_a2;

d_a1 = 0.0;
d_a2 = 0.0;

for(k=0;k<length;k++) {  

   for(n=1;n<=q_int;n++) { 
       
      d_a1 += accept_1[k][n];
      
      for(l=0;l<length;l++) {  
        
         for(n2=1;n2<=q_int;n2++) {  
             
            d_a2 += accept_2[k][l][n][n2];
            
         };
      };
   };
}; 

double d_n = ((d_a1) / (((double)(length*q_int*n1))) + (d_a2)/(double)(length*q_int*length*q_int*n1))/2.0;

std::cout << " Number of accepted steps : h : " << d_a1 << " " << (d_a1) / (((double)(length*q_int*n1))) 
<< " J : " << d_a2 << " " << (d_a2)/(double)(length*q_int*length*q_int*n1) << " Fraction : " << d_n << std::endl;

    n2 = 0;

    for(k=0;k<length;k++) {
    
	for(l=0;l<length;l++) {
	
            J_coupling[k][l] = 0.0;		

	};    
    };

    for(k=0;k<length;k++) {
        
      av_m3[k] = 0.0;   
        
      for(n=1;n<=q_int;n++) {       
        
        av_m3[k] += sqrt(pow(h[k][n],2));  
          
        for(l=0;l<length;l++) {
   
          for(n2=1;n2<=q_int;n2++) {              
          
           // if(h[k][n] != 0.0 && h[l][n2] != 0.0 && J[k][l][n][n2] != 0.0) J_coupling[k][l] += (J[k][l][n][n2]) * log((J[k][l][n][n2]) / (h[k][n]*h[l][n2]));
            
              J_coupling[k][l] += sqrt(pow(J[k][l][n][n2],2));
              
            };
          };
        };
        
       };
    
    double av_kl2=0.0,av_kl3=0.0;   
    
    for(k=0;k<length;k++) {
    
	av_k2[k] = 0.0;
        av_l2[k] = 0.0;	
    
    }

    for(k=0;k<length;k++) {
            
      for(l=0;l<length;l++) {
            
          av_k2[k]  += J_coupling[k][l]/((double)length);
          av_l2[l]  += J_coupling[k][l]/((double)length);
          av_kl2    += J_coupling[k][l]/((double)(length*length));
          
        };
       };   
    
    for(k=0;k<length;k++) {
            
      for(l=0;l<length;l++) {
      
          if(av_kl2 > 0.0) J_coupling[k][l] = J_coupling[k][l] - (av_k2[k]*av_l2[l])/av_kl2;
          
      };
    };
   
    for(k=0;k<length;k++) {
            
      for(l=0;l<length;l++) {
      
          J_coupling[k][l] = 0.5*(J_coupling[k][l] + J_coupling[l][k]);
          
      };
    };  
    
    av_kl2 = 0.0;
   
    for(k=0;k<length;k++) {

	av_k2[k] = 0.0;
        av_l2[k] = 0.0;
        av_m2[k] = 0.0;	

    }

    for(k=0;k<length;k++) {   
        
        for(l=0;l<length;l++) {
              
            for(m=0;m<length;m++) {
                
                    if(av_m3[k] > 0.0 & av_m3[l] > 0.0 && av_m3[m] > 0.0) V_coupling[k][l][m] = J_coupling[k][l]*J_coupling[k][m]*J_coupling[l][m]/(av_m3[k]*av_m3[l]*av_m3[m]);
                    
                    av_k2[k] += V_coupling[k][l][m]/((double)(length*length));
                    av_l2[l] += V_coupling[k][l][m]/((double)(length*length));
                    av_m2[m] += V_coupling[k][l][m]/((double)(length*length));
                    
                    av_kl2   += V_coupling[k][l][m]/((double)(length*length*length));

            };
           };
          };    
          
    for(k=0;k<length;k++) {   
        
        for(l=0;l<length;l++) {
            
            for(m=0;m<length;m++) {          
    
              if(av_kl2 > 0.0) V_coupling[k][l][m] = V_coupling[k][l][m] - (av_k2[k]*av_l2[l]*av_m2[m])/av_kl2;
              
            };
        };
    };
    
    for(k=0;k<length;k++) {   
        
        for(l=0;l<length;l++) {
            
            for(m=0;m<length;m++) {   
               
              V_coupling[k][l][m] = 0.33*(V_coupling[k][l][m]+V_coupling[k][m][l]+V_coupling[m][l][k]);
                
            }
        }
    }
    
    double buff_J1,buff_J2;
    int    buff_k1,buff_k2;
    int    buff_l1,buff_l2;
    int    buff_m1,buff_m2;

    int *k_index;
    int *l_index;
    double *J_index;


    int *k_index2;
    int *l_index2;
    int *m_index2;
    double *V_index;

    int r,r2;

    r = 0;
    r2 = 0;

    for(k=0;k<length;k++) {

       for(l=k;l<length;l++) {

	     r ++;

         for(m=l;m<length;m++) {

	     r2 ++;
	 }
       }
    }  

    k_index = (int *) malloc(sizeof(int)*((r+2)));
    l_index = (int *) malloc(sizeof(int)*((r+2)));
    J_index = (double *) malloc(sizeof(double)*((r+2)));

    k_index2 = (int *) malloc(sizeof(int)*((r2+2)));
    l_index2 = (int *) malloc(sizeof(int)*((r2+2)));
    m_index2 = (int *) malloc(sizeof(int)*((r2+2)));
    V_index  = (double *) malloc(sizeof(double)*((r2+2)));


    r = 0;
    r2 = 0;

    for(k=0;k<length;k++) {
        
        for(l=k;l<length;l++) {
          
            k_index[r] = k;
            l_index[r] = l;
            J_index[r] = J_coupling[k][l];
         
         //   std::cout << length << " " << k << " " << l << " " << r << " r " << std::endl;
            
            r++;
            
            for(m=l;m<length;m++) {
              
                k_index2[r2] = k;
                l_index2[r2] = l;
                m_index2[r2] = m;
                
                V_index[r2]  = V_coupling[k][l][m];
                
                r2 ++;
                
            };
        };
    };
    
   Quicksort(J_index,k_index,l_index,0,r);
   Quicksort2(V_index,k_index2,l_index2,m_index2,0,r2);    
   
    n = 1;
  
     for(k=0;k<=r;k++) {
     
       if(k_index[k] != l_index[k]) {  
         
        fprintf(fp,"%d\t%d\t%15.25e\n",k_index[k]+1,l_index[k]+1,J_index[k]);       
        
       }; 
     };
    
     for(k=0;k<=r2;k++) {     
     
       if(k_index2[k] != l_index2[k] && l_index2[k] != m_index2[k]) {   
         
          fprintf(fp2,"%d\t%d\t%d\t%15.25e\n",k_index2[k]+1,l_index2[k]+1,m_index2[k]+1,V_index[k]); 
           
       };
     };
    
   free(k_index);
   free(l_index);
   free(J_index);
   free(k_index2);
   free(l_index2);
   free(m_index2);
   free(V_index);

    fprintf(fp,"\n");
    fprintf(fp2,"\n");    
     
 };
}
    
    std::cout << " finalized " << std::endl; 
    
    fclose(fp);
    fclose(fp2);
}

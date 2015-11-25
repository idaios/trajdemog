//  Converts an arbitrary frequency trajectory into an approximating step function.
// Takes as input a trajectory in the format generated by forwardsel, and converts it to a step
//  function specified in a format appropriate for mssel.
//
//  The step function approximation depends on a series of boundary points stored in freqints.h .

#include "freqints.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main(int argc, char **argv )
{
  double a, b,  *freq , *mytimes, p ,  t , *fudge ;
  
  int popsize, ntpoints, ntrajs, gen, maxtpoints, currentinterval,  i , countints, myint, fixgeneration = -1 ;

  int  traj , numinterval  ;

  char line[1001] ;
 

  a = 1.0 ;
  b = 0.0 ;

  if( argc > 1 && strcmp( argv[1], "shift") != 0){
    if( argv[1][0] == '-' ) { 
      fprintf(stderr,"stepftn [ a b ] ,  (times are transformed to : a*t + b ) \n");
      exit(1);
    }
    else { 
      a = strtod( argv[2],NULL ) ;
      b = strtod( argv[3],NULL ) ;
    }
  }
  maxtpoints = 100000 ;
  
  freq = (double *)malloc( (unsigned)maxtpoints*sizeof(double) ) ;
  mytimes = (double *)malloc( (unsigned)maxtpoints*sizeof(double) ) ;
  
  do {
    if( fgets( line, 1000, stdin) == NULL ) {
      fprintf(stderr,"stepftn [ a b ] ,  (times are transformed to : a*t + b ) \n");
      exit(1);
    }
    if( line[0] == '/' ) printf("%s",line);
  }while ( line[0] == '/' );
 

  sscanf(line," %d",&ntrajs);

  printf("ntraj: %d\n", ntrajs);

  printf("npop: 1\n");

  fgets( line, 1000, stdin);

  if( line[0] != '#' ) { fprintf(stderr," file format error.\n"); exit(0); }

  

  for( traj = 0; traj < ntrajs; traj++){
  
    countints = 0;

    fixgeneration = -1;
  
    while( 1 ) {
      
      if( fgets( line, 1000, stdin) == NULL ){
	break;
      }
      
      if( line[0] == '#' ) break; 
   
      sscanf(line,"%lf %lf",&t,  &p);

      if( p > 0.99999999 && fixgeneration == -1 && argc > 1 && strcmp(argv[1], "shift")== 0) 
	fixgeneration = countints;
      
      if( countints > maxtpoints -2 ){
	maxtpoints += 5000 ;
	freq = (double *)realloc(freq,(unsigned)(maxtpoints*sizeof(double)) );
	mytimes = (double *)realloc(mytimes,(unsigned)(maxtpoints*sizeof(double)) );	 
      }
    
  
	
	
      if( countints == 0 ||  mytimes[countints-1] != t ){
	mytimes[countints] = t;
	freq[countints]  =  p ;
	countints++;
      }
      
    }
    if( argc > 1 && strcmp(argv[1], "shift")== 0 ) 
    	printf("n:\t%d\n", fixgeneration+1); //countints-1);
    else
	printf("n:\t%d\n", countints-1);
    
    for( i= countints-1; i>=0 ; i--) 
      {
	
	if( argc > 1 && strcmp(argv[1], "shift")== 0 )
	{
		if( fixgeneration + i - countints+1 >= 0)
	  		printf("%lf\t%lf\n", mytimes[countints-1] - mytimes[i], freq[fixgeneration + i - countints+1] ) ; 
	}
	else
	{
	  printf("%lf\t%lf\n", a*( mytimes[countints-1]-mytimes[i]) + b , freq[i-1]);
	}	
	
      }
  }
  
}
 

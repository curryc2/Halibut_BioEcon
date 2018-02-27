DATA_SECTION
 init_int Fyear;
 init_int Lyear;
 init_vector Catch(Fyear,Lyear);
 init_vector CPUE(Fyear,Lyear+1);
 !! cout << Catch << endl;                                       // Check that the data are read in
 !! cout << CPUE << endl;                                       // Check that the data are read in

PARAMETER_SECTION
 init_number logr;
 init_number logq;
 init_number logK;
 number r;
 number q;
 number K;
 vector Bio(Fyear,Lyear+1);
 
 objective_function_value objn;

INITIALIZATION_SECTION
 logr -0.6
 logq -9
 logK 8.5 

PROCEDURE_SECTION
 int Year;
 dvariable SS;
 
 // Convert from log to normal space
 r = mfexp(logr);
 q = mfexp(logq);
 K = mfexp(logK);

 // Project the model forward
 Bio(Fyear) = K;
 for (Year=Fyear;Year<=Lyear;Year++)
  {
   Bio(Year+1) = Bio(Year) + r*Bio(Year)*(1.0-Bio(Year)/K) - Catch(Year);
   if (Bio(Year+1) < 1) Bio(Year+1) = 1;
  } 
  
 // Compute the sum of squares  
 SS = 0;
 for (Year=Fyear;Year<=Lyear+1;Year++)
  if (CPUE(Year) > 0)
   SS += square(log(q*Bio(Year)/CPUE(Year)));
 objn = SS;  
 cout << SS << endl;
 
REPORT_SECTION
 report << mfexp(logr) << endl;
 report << Bio(Lyear+1) << endl;
 report << " " << endl;
 report << Bio << endl;
 

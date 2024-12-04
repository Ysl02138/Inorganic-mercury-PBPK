Integrate (Lsodes, 1e-8, 1e-8 , 1);

Simulation { #single oral dose of 1000 ug Hg/kg/d
  
  BW0 = 0.02;
  BWgrowth= 0;
  sex= 2;
  Growthrate = 0.0;
  TChng 	= 0.05;

  PDose = PerDose(1000.0, 24,  0, 0.05);
  IVDose = PerDose(0.0, 24,  0, 0.003);
  Drink = 0;
  
  expowk = PerDose(1.0, 168,  0, 0.05);
  expodur = PerDose(1.0, 1344, 0, 0.05);
  
  PrintStep (BalanceCheck, KGavage, CBldU, CBrnU, 0, 480, 1);

 } 
  
END.




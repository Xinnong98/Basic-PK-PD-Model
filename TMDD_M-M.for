**********************************************************************
C                           ADAPT                                     *
C                         Version 5                                   *
C**********************************************************************
C                                                                     *
C                           MODEL                                     *
C                                                                     *
C    This file contains Fortran subroutines into which the user       *
C    must enter the relevant model equations and constants.           *
C    Consult the User's Guide for details concerning the format for   *
C    entered equations and definition of symbols.                     *
C                                                                     *
C       1. Symbol- Parameter symbols and model constants              *
C       2. DiffEq- System differential equations                      *
C       3. Output- System output equations                            *
C       4. Varmod- Error variance model equations                     *
C       5. Covmod- Covariate model equations (population)             *
C       6. Prior - Parameter mean and covariance values               *
C       7. Sparam- Secondary parameters                               *
C       8. Amat  - System state matrix                                *
C                                                                     *
C**********************************************************************

C######################################################################C

        Subroutine SYMBOL
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

CC
C----------------------------------------------------------------------C
C   Enter as Indicated                                                 C
C----c-----------------------------------------------------------------C

      NDEqs   = 3   ! Enter # of Diff. Eqs.
      NSParam = 6   ! Enter # of System Parameters.
      NVparam = 2   ! Enter # of Variance Parameters.
      NSecPar = 0   ! Enter # of Secondary Parameters.
      NSecOut = 0   ! Enter # of Secondary Outputs (not used).
      Ieqsol  = 1   ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' MM '

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C
     
	   Psym(1)='kel'
	   Psym(2)='kt'
	   Psym(3)='kon'
	   Psym(4)='kint'
	   Psym(5)='Vc'
	   Psym(6)='kdeg'
      
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}    C
C----c-----------------------------------------------------------------C
      
       PVsym(1)='sigma1'
       PVsym(2)='sigma2'
     
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Secondary Parameter {eg: PSsym(1)='CLt'}     C
C----c-----------------------------------------------------------------C



C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End
        
C######################################################################C

      
        Subroutine DIFFEQ(T,X,XP)
        Implicit None
 
        Include 'globals.inc'
        Include 'model.inc'

        Real*8 T,X(MaxNDE),XP(MaxNDE)
        Real*8 kel,kt,kon,koff,kint,Vc,kdeg,ksyn,Km,Vm

     
        
CC   
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C
      
       kel=P(1) !1/h
	   kt=P(2)  !1/h
	   kon=P(3) !1/nM/h secong order
	   koff=0.1*kon !1/h where 0.1 is KD (nM)
	   kint=P(4)!1/h
	   Vc=P(5)  !L/kg
	   kdeg=P(6)!1/h
	   ksyn=kdeg*IC(3) !nmol/h where IC(3) is amount
	   Km = (kint+koff)/kon !nM
	   Vm = kint*X(3) !nmol/h
	   
	   XP(1)=-(kel+kt)*X(1)-(Vm*X(1))/((Km*Vc)+X(1))+kt*X(2)	!Free Drug Amount in Central Compartment (nmol)   
	   XP(2)=-kt*X(2)+kt*X(1) !Drug Amount distributed into Tissue Compartment (nmol)
	   XP(3)=ksyn-kdeg*X(3)-(kint-kdeg)*(X(3)*X(1)/((Km*Vc)+X(1))) !Total Receptor Amount (nmol)
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End



C######################################################################C

        Subroutine OUTPUT(Y,T,X)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 Y(MaxNOE),T,X(MaxNDE)
        Real*8 kel,kt,kon,koff,kint,Vc,kdeg,ksyn,Km,Vm

        
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C
       kel=P(1)  !1/h
	   kt=P(2)   !1/h
	   kon=P(3)  !1/nM/h secong order
	   koff=0.1*kon !1/h
	   kint=P(4) !1/h
	   Vc=P(5)   !L/kg
	   kdeg=P(6) !1/h
	   ksyn=kdeg*IC(3) !nmol/h zero order IC(3) is amount
	   Km = (kint+koff)/kon !nM
	   Vm = kint*X(3)  !nmol/h

	   

       Y(1)=X(1)/Vc*19.7 !Output is Concentration (where X(1) is Amount
      
      
C----------------------------------------------------------------------C
C---------------------------------------------------------------------C
C
        Return
        End


C######################################################################C

        Subroutine VARMOD(V,T,X,Y)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 V(MaxNOE),T,X(MaxNDE),Y(MaxNOE)

CC
C----------------------------------------------------------------------C
C   Enter Variance Model Equations Below                               C
C         {e.g. V(1) = (PV(1) + PV(2)*Y(1))**2 }                       C
C----c-----------------------------------------------------------------C
  
      V(1) = (PV(1) + PV(2)*Y(1))**2 
  
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End


C######################################################################C

        Subroutine PRIOR(Pmean,Pcov,ICmean,ICcov)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)
        Real*8 Pcov(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcov(MaxNDE,MaxNDE)

CC
C----------------------------------------------------------------------C
C  Enter Nonzero Elements of Prior Mean Vector                         C
C          {  e.g. Pmean(2) = 10.0    }                                C
C----c-----------------------------------------------------------------C



C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
CC
C----------------------------------------------------------------------C
C   Enter Nonzero Elements of Covariance Matrix (Lower Triang.)       C
C         {  e.g. Pcov(2,1) = 0.25    }                                C
C----c-----------------------------------------------------------------C



C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End
        
C######################################################################C
	  
        Subroutine SPARAM(PS,P,IC)
        Implicit None

        Include 'globals.inc'

        Real*8 PS(MaxNSECP), P(MaxNSP+MaxNDE), IC(MaxNDE) 

CC
C----------------------------------------------------------------------C
C   Enter Equations Defining Secondary Paramters                       C
C           {  e.g.  PS(1) = P(1)*P(2)   }                             C
C----c-----------------------------------------------------------------C

        
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End


C######################################################################C
	  

        Subroutine AMAT(A)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 A(MaxNDE,MaxNDE)

        DO I=1,Ndeqs
           Do J=1,Ndeqs
              A(I,J)=0.0D0
           End Do
        End Do

CC
C----------------------------------------------------------------------C
C   Enter non zero elements of state matrix  {e.g.  A(1,1) = -P(1) }   C
C----c-----------------------------------------------------------------C



C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C
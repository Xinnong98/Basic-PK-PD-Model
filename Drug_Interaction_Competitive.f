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
C       1. Symbol-  Parameter symbols and model constants             *
C       2. DiffEq-  System differential equations                     *
C       3. Output-  System output equations                           *
C       4. Varmod-  Error variance model equations                    *
C       5. Covmod-  Covariate model equations (ITS,MLEM)              *
C       6. Popinit- Population parameter initial values (ITS,MLEM)    *
C       7. Prior -  Parameter mean and covariance values (ID,NPD,STS) *
C       8. Sparam-  Secondary parameters                              *
C       9. Amat  -  System state matrix                               *
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

      NDEqs   =  14   ! Enter # of Diff. Eqs.
      NSParam =  16   ! Enter # of System Parameters.
      NVparam =  2  ! Enter # of Variance Parameters.
      NSecPar =  0  ! Enter # of Secondary Parameters.
      NSecOut =  0  ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1  ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' interaction '

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C

	   Psym(1)='kinss'
	   Psym(2)='kout'
	   Psym(3)='Sgmax'
	   Psym(4)='SC50g'
	   Psym(5)='Sbmax'
       Psym(6)='SC50b'
	   Psym(7)='kd'
	   Psym(8)='ga1'
	   Psym(9)='ga2'
	   Psym(10)='psi'
	   Psym(11)='Igmax'
	   Psym(12)='IC50g'
	   Psym(13)='Ibmax'
       Psym(14)='IC50b'
	   Psym(15)='psi2'
	   Psym(16)='KT'
	   



CC
C----------------------------------------------------------------------C
C   Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}    C
C----c-----------------------------------------------------------------C
       PVsym(1)='int'
	   PVsym(2)='sig'


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
		Real*8 kout,Sgmax,SC50g,Sbmax,SC50b,kd,kinss,ga1,ga2,psi
		Real*8 gf1,gf2,gf3,bf1,bf2,bf3,Igmax,IC50g,Ibmax,IC50b
		Real*8 GF4,GF5,GF6,BF4,BF5,BF6,psi2,KT
CC
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C
		 kinss=P(1)
		 kout=P(2)
		 Sgmax=P(3)
		 SC50g=P(4)
		 Sbmax=P(5)
		 SC50b=P(6)
		 kd=P(7)
		 ga1=P(8)
		 ga2=P(9)
		 psi=P(10)
		 Igmax=P(11)
		 IC50g=P(12)
		 Ibmax=P(13)
		 IC50b=P(14)
		 psi2=p(15)
		 KT=P(16)
		 !COMPETITIVE
		 gf1 = (6/(psi*SC50g))**ga1
		 gf2 = (10/(psi*SC50g))**ga1
		 gf3 = (20/(psi*SC50g))**ga1
		 gf4 = (6/(psi2*IC50g))**ga1
		 gf5 = (10/(psi2*IC50g))**ga1
		 gf6 = (20/(psi2*IC50g))**ga1
		 
		 bf1 = (50/(psi*SC50b))**ga2
		 bf2 = (200/(psi*SC50b))**ga2
		 bf3 = (500/(psi*SC50b))**ga2
		 bf4 = (50/(psi2*IC50b))**ga2
		 bf5 = (200/(psi2*IC50b))**ga2
		 bf6 = (500/(psi2*IC50b))**ga2

		 XP(1)=X(5)-kout*X(1)
		 XP(2)=X(5)*(1-(Igmax*gf4+Ibmax*bf4)/(bf4+gf4+1))-
     C      kout*(1+X(8))*X(2)
		 XP(3)=X(5)*(1-(Igmax*gf5+Ibmax*bf5)/(bf5+gf5+1))-
     C      kout*(1+X(11))*X(3)
		 XP(4)=X(5)*(1-(Igmax*gf6+Ibmax*bf6)/(bf6+gf6+1))-
     C      kout*(1+X(14))*X(4)
		
		 XP(5)=kd*X(5)*(1-X(5)/kinss)
		 
		 XP(6)=kt*((Sgmax*gf1+Sbmax*bf1)/(bf1+gf1+1)-X(6))
		 XP(7)=kt*(X(6)-X(7))
		 XP(8)=kt*(X(7)-X(8))
		 XP(9)=kt*((Sgmax*gf2+Sbmax*bf2)/(bf2+gf2+1)-X(9))
		 XP(10)=kt*(X(9)-X(10))
		 XP(11)=kt*(X(10)-X(11))
		 XP(12)=kt*((Sgmax*gf3+Sbmax*bf3)/(bf3+gf3+1)-X(12))
		 XP(13)=kt*(X(12)-X(13))
		 XP(14)=kt*(X(13)-X(14))
		 
		 
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

        Real*8 Y(MaxNOE),T,X(MaxNDE),V


CC
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C
		 Y(1)=X(1)
		 Y(2)=X(2)
		 Y(3)=X(3)
		 Y(4)=X(4)

		 
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
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
		 V(2) = (PV(1) + PV(2)*Y(2))**2
		 V(3) = (PV(1) + PV(2)*Y(3))**2
		 V(4) = (PV(1) + PV(2)*Y(4))**2

		 
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine COVMOD(Pmean, ICmean, PC)
C  Defines any covariate model equations (MLEM, ITS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 PC(MaxNCP)
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)

CC
C----------------------------------------------------------------------C
C     Enter # of Covariate Parameters                                  C
C----c-----------------------------------------------------------------C

        NCparam = 0    ! Enter # of Covariate Parameters.

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Covariate Params {eg: PCsym(1)='CLRenal'}         C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   For the Model Params. that Depend on Covariates Enter the Equation C
C         {e.g. Pmean(1) =  PC(1)*R(2) }                               C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine POPINIT(PmeanI,ICmeanI,PcovI,ICcovI, PCI)
C  Initial parameter values for population program parameters (ITS, MLEM)

        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 PmeanI(MaxNSP+MaxNDE), ICmeanI(MaxNDE)
        Real*8 PcovI(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcovI(MaxNDE,MaxNDE)
        Real*8 PCI(MaxNCP)

CC
C----------------------------------------------------------------------C
C  Enter Initial Values for Population Means                           C
C          {  e.g. PmeanI(1) = 10.0    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Initial Values for Pop. Covariance Matrix (Lower Triang.)    C
C         {  e.g. PcovI(2,1) = 0.25    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Values for Covariate Model Parameters                        C
C         {  e.g. PCI(1) = 2.0    }                                    C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End

C######################################################################C

        Subroutine PRIOR(Pmean,Pcov,ICmean,ICcov)
C  Parameter mean and covariance values for MAP estimation (ID,NPD,STS)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Integer I,J
        Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)
        Real*8 Pcov(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcov(MaxNDE,MaxNDE)

CC
C----------------------------------------------------------------------C
C  Enter Nonzero Elements of Prior Mean Vector                         C
C          {  e.g. Pmean(1) = 10.0    }                                C
C----c-----------------------------------------------------------------C


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
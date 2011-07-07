      subroutine x_rks_s(r,f,dfdra)
      implicit none
c
c     This subroutine evaluates the spin polarised exchange functional
c     in the Local Density Approximation [1], and the corresponding
c     potential. Often this functional is referred to as the Dirac 
c     functional [2] or Slater functional.
c
c     [1] F. Bloch, Zeitschrift fuer Physik, Vol. 57 (1929) 545.
c
c     [2] P.A.M. Dirac, Proceedings of the Cambridge Philosophical
c         Society, Vol. 26 (1930) 376.
c
c     Parameters:
c
c     r     the total electron density
c     f     On return the functional value
c     dfdra On return the derivative of f with respect to alpha electron
c           density
c
c
c     Ax = -3/4*(6/pi)**(1/3)
c     Bx = -(6/pi)**(1/3)
c     C  = (1/2)**(1/3)
      real*8 Ax, Bx, C
      parameter(Ax    = -0.930525736349100025002010218071667d0)
      parameter(Bx    = -1.24070098179880003333601362409556d0)
      parameter(C     =  0.793700525984099737375852819636154d0)
c
      real*8 third
      parameter(third =  0.333333333333333333333333333333333d0)
c
      real*8 r
      real*8 f, dfdra
c
      real*8 ra13
c
      ra13  = C*r**third
      f     = Ax*r*ra13
      dfdra = Bx*ra13 
c
      end
c-----------------------------------------------------------------------
      subroutine x_uks_s(ra,rb,f,dfdra,dfdrb)
      implicit none
c
c     This subroutine evaluates the spin polarised exchange functional
c     in the Local Density Approximation [1], and the corresponding
c     potential. Often this functional is referred to as the Dirac 
c     functional [2] or Slater functional.
c
c     [1] F. Bloch, Zeitschrift fuer Physik, Vol. 57 (1929) 545.
c
c     [2] P.A.M. Dirac, Proceedings of the Cambridge Philosophical
c         Society, Vol. 26 (1930) 376.
c
c     Parameters:
c
c     ra    the alpha electron density
c     rb    the beta  electron density
c     f     On return the functional value
c     dfdra On return the derivative of f with respect to ra
c     dfdrb On return the derivative of f with respect to rb
c
c
c     Ax = -3/4*(6/pi)**(1/3)
c     Bx = -(6/pi)**(1/3)
      real*8 Ax, Bx
      parameter(Ax    = -0.930525736349100025002010218071667d0)
      parameter(Bx    = -1.24070098179880003333601362409556d0)
c
      real*8 third
      parameter(third =  0.333333333333333333333333333333333d0)
c
      real*8 ra, rb
      real*8 f, dfdra, dfdrb
c
      real*8 ra13, rb13
c
      ra13  = ra**third
      rb13  = rb**third
      f     = Ax*(ra*ra13+rb*rb13) 
      dfdra = Bx*ra13 
      dfdrb = Bx*rb13 
c
      end
      subroutine c_rks_vwn5(r,f,dfdra)
      implicit none
c
c     This subroutine evaluates the Vosko, Wilk and Nusair correlation
c     functional number 5 [1] for the closed shell case, with the 
c     parametrisation as given in table 5.
c
c     The original code was obtained from Dr. Phillip Young,
c     with corrections from Dr. Paul Sherwood.
c
c     [1] S.H. Vosko, L. Wilk, and M. Nusair
c         "Accurate spin-dependent electron liquid correlation energies
c          for local spin density calculations: a critical analysis",
c         Can.J.Phys, Vol. 58 (1980) 1200-1211.
c
c     Parameters:
c
c     r      the total electron density
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to the alpha
c            electron density
c
      real*8 r
      real*8 f, dfdra
c
      real*8 t4,t5,t6,t7
      real*8 a2,b2,c2,d2
      real*8 P1,P2,P3,P4
      real*8 alpha_rho13
      real*8 srho,srho13
      real*8 iv,iv2,inv,i1,i2,i3
      real*8 pp1,pp2
c
      real*8 n13, n43, n16, one, two, four
      parameter(n13 = 0.333333333333333333333333333333d0)
      parameter(n43 = 1.333333333333333333333333333333d0)
      parameter(n16 = 0.166666666666666666666666666666d0)
      parameter(one = 1.0d0)
      parameter(two = 2.0d0)
      parameter(four= 4.0d0)
C 
C VWN interpolation parameters
C
C paramagnetic
      a2 =  0.0621814d0 
      b2 =  3.72744d0
      c2 = 12.9352d0    
      d2 = -0.10498d0
C
C t4 = (1/(4/3)*pi)**(1/3)
      t4=0.620350490899399531d0
C
C t5 = 0.5/(2**(1/3)-1)
      t5 = 1.92366105093153617d0
C
C t6 = 2/(3*(2**(1/3)-1))
      t6 = 2.56488140124204822d0
C
C t7 = 2.25*(2**(1/3)-1)
      t7 = 0.584822362263464735d0
C
C Paramagnetic interpolation constants
C 
      P1 = 6.15199081975907980d0
      P2 = a2 *0.5d0
      P3 = 0.193804554230887489d-02 *0.5d0
      P4 = 0.775665897562260176d-01 *0.5d0
C
C closed shell case
      srho        = r
      srho13      = srho**n13
      alpha_rho13 = (0.5d0**n13)*srho
      iv2         = T4/srho13
      iv          = sqrt(iv2)
C
C paramagnetic
      inv = 1.0d0/(iv2+b2*iv+c2)
      i1  = log(iv2*inv)
      i2  = log((iv-d2)*(iv-d2)*inv)
c corrected b1->b2 ps Apr98
      i3  = atan(P1/(2.0d0*iv+b2))
      pp1 = P2*i1 + P3*i2 + P4*i3
      pp2 = a2*(1.0d0/iv-iv*inv*(1.0d0+b2/(iv-d2))) 
C
      f     = pp1*srho
      dfdra = pp1 - n16*iv*pp2
c
      end
c-----------------------------------------------------------------------
      subroutine c_uks_vwn5(ra,rb,f,dfdra,dfdrb)
      implicit none
c
c     This subroutine evaluates the Vosko, Wilk and Nusair correlation
c     functional number 5 [1], with the parametrisation as given in
c     table 5.
c
c     The original code was obtained from Dr. Phillip Young,
c     with corrections from Dr. Paul Sherwood.
c
c     [1] S.H. Vosko, L. Wilk, and M. Nusair
c         "Accurate spin-dependent electron liquid correlation energies
c          for local spin density calculations: a critical analysis",
c         Can.J.Phys, Vol. 58 (1980) 1200-1211.
c
c     Parameters:
c
c     ra     the alpha-electron density
c     rb     the beta-electron density
c     f      On return the functional value
c     dfdra  On return the derivative of f with respect to ra
c     dfdrb  On return the derivative of f with respect to rb
c
      real*8 ra, rb
      real*8 f, dfdra, dfdrb
c
      real*8 t1,t2,t4,t5,t6,t7
      real*8 a1,b1,c1,d1
      real*8 a2,b2,c2,d2
      real*8 a3,b3,c3,d3
      real*8 S1,S2,S3,S4
      real*8 P1,P2,P3,P4
      real*8 F1,F2,F3,F4
      real*8 inter1,inter2
      real*8 alpha_rho13,beta_rho13
      real*8 srho,srho13
      real*8 iv,iv2,inv,i1,i2,i3
      real*8 vwn1,vwn2
      real*8 ss1,ss2,pp1,pp2,ff1,ff2
      real*8 zeta,zeta3,zeta4
      real*8 tau,dtau,v
c
      real*8 n13, n43, n16, one, two, four, tn13, tn43
c     tn13 = 2**(1/3)
c     tn43 = 2**(4/3)
      parameter(n13 = 0.333333333333333333333333333333d0)
      parameter(n43 = 1.333333333333333333333333333333d0)
      parameter(n16 = 0.166666666666666666666666666667d0)
      parameter(tn13= 1.25992104989487316476721060727823d0)
      parameter(tn43= 2.51984209978974632953442121455646d0)
      parameter(one = 1.0d0)
      parameter(two = 2.0d0)
      parameter(four= 4.0d0)
C 
C VWN interpolation parameters
C
C spin stiffness
      a1 = -0.337737278807791058d-01
      b1 =  1.13107d0
      c1 = 13.0045d0
      d1 = -0.00475840d0
C paramagnetic
      a2 =  0.0621814d0 
      b2 =  3.72744d0
      c2 = 12.9352d0    
      d2 = -0.10498d0
C ferromagnetic
c try cadpac/nwchem value (.5*a2)
      a3 = 0.0310907d0
      b3 =  7.06042d0
      c3 = 18.0578d0
      d3 = -0.32500d0
C
C t4 = (1/(4/3)*pi)**(1/3)
      t4=0.620350490899399531d0
C
C t5 = 0.5/(2**(1/3)-1)
      t5 = 1.92366105093153617d0
C
C t6 = 2/(3*(2**(1/3)-1))
      t6 = 2.56488140124204822d0
C
C t7 = 2.25*(2**(1/3)-1)
      t7 = 0.584822362263464735d0
C
C Spin stiffness interpolation constants
C
      S1 = 7.12310891781811772d0
      S2 = a1 *0.5d0
      S3 = -0.139834647015288626d-04 *0.5d0
      S4 = -0.107301836977671539d-01 *0.5d0
C
C Paramagnetic interpolation constants
C 
      P1 = 6.15199081975907980d0
      P2 = a2 *0.5d0
      P3 = 0.193804554230887489d-02 *0.5d0
      P4 = 0.775665897562260176d-01 *0.5d0
C
C Ferromagnetic interpolation constants
C
      F1 = 4.73092690956011364d0
      F2 = a3 *0.5d0
c
c      F3 = -0.244185082989490298d-02 *0.5d0
c      F4 = -0.570212323620622186d-01 *0.5d0
c
c  try nwchem values
c
      F3 = 0.224786709554261133d-02 
      F4 = 0.524913931697809227d-01 
C
C Interpolation intervals
C 
      inter1 =  1.0d0-1.0d-10
      inter2 = -1.0d0+1.0d-10
C
C open shell case
      alpha_rho13 = ra**n13
      beta_rho13  = rb**n13
      srho        = ra+rb
      srho13      = srho**n13
      iv2         = T4/srho13
      iv          = sqrt(iv2)
C
C spin-stiffness
      inv = 1.0d0/(iv2+b1*iv+c1)
      i1  = log(iv2*inv)
      i2  = log((iv-d1)*(iv-d1)*inv)
      i3  = atan(S1/(2.0d0*iv+b1))
      ss1 = S2*i1 + S3*i2 + S4*i3
      ss2 = a1*(1.0d0/iv-iv*inv*(1.0d0+b1/(iv-d1)))
C
C paramagnetic
      inv = 1.0d0/(iv2+b2*iv+c2)
      i1  = log(iv2*inv)
      i2  = log((iv-d2)*(iv-d2)*inv)
c corrected b1->b2 ps Apr98
      i3  = atan(P1/(2.0d0*iv+b2))
      pp1 = P2*i1 + P3*i2 + P4*i3
      pp2 = a2*(1.0d0/iv-iv*inv*(1.0d0+b2/(iv-d2))) 
C
C ferromagnetic
      inv = 1.0d0/(iv2+b3*iv+c3)
      i1  = log(iv2*inv)
      i2  = log((iv-d3)*(iv-d3)*inv)
      i3  = atan(F1/(2.0d0*iv+b3))
      ff1 = F2*i1 + F3*i2 + F4*i3
      ff2 = a3*(1.0d0/iv-iv*inv*(1.0d0+b3/(iv-d3)))
C
C polarisation function
c
      zeta  = (ra-rb)/srho
      zeta3 = zeta*zeta*zeta
      zeta4 = zeta3*zeta
      if(zeta.gt.inter1) then
         vwn1 = (tn43-two)*t5
         vwn2 = (tn13)*t6
      elseif(zeta.lt.inter2) then
         vwn1 = (tn43-two)*t5
         vwn2 = -(tn13)*t6
      else
         vwn1 = ((one+zeta)**n43+(one-zeta)**n43-two)*t5
         vwn2 = ((one+zeta)**n13-(one-zeta)**n13)*t6
      endif
      ss1  = ss1*t7
      ss2  = ss2*t7 
      tau  = ff1-pp1-ss1
      dtau = ff2-pp2-ss2
c
      v = pp1+vwn1*(ss1+tau*zeta4)
      f = v*srho
c
      t1 = v - n16*iv*(pp2+vwn1*(ss2+dtau*zeta4))
      t2 = vwn2*(ss1+tau*zeta4)+vwn1*four*tau*zeta3
      dfdra = t1+t2*(one-zeta)
      dfdrb = t1-t2*(one+zeta)
c
      end

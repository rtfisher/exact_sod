      program exact_sod 

*/////////////////////////////////////////////////////////////////
*/
*/  Program calculates the exact solution to Sod-shock tube class
*/  problems -- namely shock tubes which produce shocks, contact
*/  discontinuities, and rarefraction waves.
*/
*/  Solution is computed at locations x at time t. (Though   
*/  due to self-similarity, the exact solution is identical for
*/  identical values of x/t). Output to the file 'exact_sod.out'
*/  is in the format
*/
*/       x position, density, velocity, pressure
*/
*/  NOTE : Since the post-shock flow is nonadiabatic, whereas
*/  the flow inside the rarefaction fan is adiabatic, the problem
*/  is not left-right symmetric. In particular, the high-density
*/  initial state MUST BE input on the left side.
*/ 
*/  Written by Robert Fisher, 12/5/96.
*/
*/  Edited : 12/15/96, typo in output of pressure corrected. R.F.
*/            6/08/22, removed 'pause' command, added output comment. R.F.
*/////////////////////////////////////////////////////////////////

      implicit none

      real*8 gamma, mu2, x, t, xmax
      integer numcells

C xmax determines the size of the computational domain (-xmax, +xmax).
C numcells determines the number of cells in the output table.          

      parameter (gamma   =    1.4D0)
      parameter (mu2     =   (gamma - 1.D0) / (gamma + 1.D0))
      parameter (xmax 	 =     1.D0)
      parameter (numcells =	500)
     
      real*8 pl, pr, rhol, rhor, cl, cr, pm, pressure,
     &        rhoml, vs, vt, rhomr, vm, density, velocity 
      integer i
 
      common/block1/ pl, pr, rhol, rhor, cl, cr

      real*8 rtbis          

      External rtbis

* Define the time of the problem.

      t = 2.45D-1
 
* Define the Sod problem initial conditions for the left and right states.

      pl = 1.D0
      pr = 1.D-1

      rhol = 1.D0
      rhor = 1.25D-1      

* Define sound speeds for the left and right sides of tube.

      cl = dsqrt (gamma * pl / rhol)
      cr = dsqrt (gamma * pr / rhor)

* Solve for the postshock pressure pm.

      pm = rtbis (pr, pl, 1.D-16)

* Define the density to the left of the contact discontinuity rhoml.
 
      rhoml = rhol * (pm / pl) ** (1.D0 / gamma)

* Define the postshock fluid velocity vm.

      vm = 2.D0 * cl / (gamma - 1.D0) * (1.D0 - (pm / pl) **
     &     ( (gamma - 1.D0) / (2.D0 * gamma) )) 

* Define the postshock density rhomr.

      rhomr = rhor *  ( (pm + mu2 * pr) / (pr + mu2 * pm) )

* Define the shock velocity vs.

      vs = vm / (1.D0 - rhor / rhomr) 

* Define the velocity of the rarefraction tail, vt.

      vt = cl - vm / (1.D0 - mu2) 

* Output tables of density, velocity, and pressure at time t.

      open (unit = 6, file = 'exact_sod.out')

      do i = 0, numcells
   
        x = - xmax +  2.D0 * xmax * i / numcells
 
        if (x .le. - cl * t) then
          density = rhol
        else if (x .le. -vt * t) then
          density = rhol * (-mu2 * (x / (cl * t) ) + (1 - mu2) ) **
     &    (2.D0 / (gamma - 1.D0)) 
        else if (x .le. vm * t) then
          density = rhoml
        else if (x .le. vs * t) then
          density = rhomr
        else
          density = rhor
        end if

        if (x .le. - cl * t) then
          pressure = pl
        else if (x .le. -vt * t) then
          pressure = pl * (-mu2 * (x / (cl * t) ) + (1 - mu2) ) **
     &    (2.D0 * gamma / (gamma - 1.D0))
        else if (x .le. vs * t) then
          pressure = pm 
        else            
          pressure = pr
        end if

        if (x .le. -cl * t) then
          velocity = 0.0
        else if (x .le. -vt * t) then
          velocity = (1 - mu2) * (x / t + cl)
        else if (x .le. vs * t) then
          velocity = vm
        else 
          velocity = 0.0
        end if

        write (6, 10) x, density, velocity, pressure

      end do 

      close (6)

  10  format (E22.16, ' ', E22.16, ' ', E22.16, ' ', E22.16)  
                    
      End



      function func (pm)

*//////////////////////////////////////////////////////////////////////
*/
*/  func is obtained from an identity matching the post-shocked      
*/  pressure to the post-rarefraction pressure (true since there is
*/  no pressure jump across the contact discontinuity). We use it to
*/  numerically solve for pm given the left and right initial states.
*/
*//////////////////////////////////////////////////////////////////////
 

      implicit none

      real*8 func, pm 
      real*8 gamma, mu2
 
      parameter (gamma   =    1.4D0)
      parameter (mu2     =   (gamma - 1.D0) / (gamma + 1.D0))
   
      real*8 pl, pr, rhol, rhor, cl, cr
 
      common/block1/ pl, pr, rhol, rhor, cl, cr


      func = -2*cl*(1 - (pm/pl)**((-1 + gamma)/(2*gamma)))/
     &    -   (cr*(-1 + gamma)) + 
     &    -  (-1 + pm/pr)*((1 - mu2)/(gamma*(mu2 + pm/pr)))**0.5
      return

      end 
  
      FUNCTION rtbis(x1,x2,xacc)

*/////////////////////////////////////////////////////////////////////////
*/
*/ rtbis is borrowed from Numerical Recipes. It is a bisection algorithm,
*/ which we use to solve for pm using a call to func.
*/
*/ Note that the arguments to rtbis have been altered and the value of
*/ JMAX increased. Otherwise, it is identical to the NR version.
*/
*/////////////////////////////////////////////////////////////////////////

      INTEGER JMAX
      REAL*8 rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=100)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) then
        print *, 'root must be bracketed in rtbis'
        stop
      endif
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*5.D-1
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      print *, 'too many bisections in rtbis'
      END

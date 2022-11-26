      program kinhsh2
      implicit double precision (a-h,o-z)
      dimension aux(2),ak(2,4)
      fun1(ti,u)=u
      fun2(ti,u)=-(12.774862738556573d0*1.5d0*0.1d0*0.06d0*1.2d0*u)
     &/1000.d0

      n=0
      deltat=1.d0
c

      V= 14
      x=0
      OPEN(2,FILE='tp.txt')     !
      write(2,*)n,x
c
      time=1.d0
      do ktimestep=1,10

c-step 1:
         ak(1,1) = deltat*fun1(time,V)
         ak(2,1) = deltat*fun2(time, V)
         aux(1) = time+0.5d0* deltat
         aux(2) = V+0.5d0*ak(2,1)
c-step 2:
         ak(1,2) = deltat*fun1(aux(1),aux(2))
         ak(2,2) = deltat*fun2(aux(1),aux(2))
         aux(1) = time+0.5d0* deltat
         aux(2) = V+0.5d0*ak(2,2)
c-step 3:
         ak(1,3) = deltat*fun1(aux(1),aux(2))
         ak(2,3) = deltat*fun2(aux(1),aux(2))
         aux(1) = time+ deltat
         aux(2) = V+ak(2,3)
c-step 4:
         ak(1,4) = deltat*fun1(aux(1),aux(2))
         ak(2,4) = deltat*fun2(aux(1),aux(2))
c-synthesis:
      x=x+(ak(1,1)+2.d0*ak(1,2)+2.d0*ak(1,3)+ak(1,4))/6.d0
      v=v+(ak(2,1)+2.d0*ak(2,2)+2.d0*ak(2,3)+ak(2,4))/6.d0

      write(2,*)ktimestep,x
      time=time+deltat
      enddo
c
      CLOSE(2)
c
      end

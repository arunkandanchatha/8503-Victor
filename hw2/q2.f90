
MODULE nr
    use nrtype
    implicit none

contains
    FUNCTION brent(ax,bx,cx,func,tol,xmin)
        USE nrtype; USE nrutil, ONLY : nrerror
        IMPLICIT NONE
        REAL(dp), INTENT(IN) :: ax,bx,cx,tol
        REAL(dp), INTENT(OUT) :: xmin
        REAL(dp) :: brent
        PROCEDURE(template_function), POINTER :: func
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(dp), PARAMETER :: CGOLD=0.3819660_dp,ZEPS=1.0e-3_dp*epsilon(ax)
        INTEGER(I4B) :: iter
        REAL(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
        a=min(ax,cx)
        b=max(ax,cx)
        v=bx
        w=v
        x=v
        e=0.0
        fx=func(x)
        fv=fx
        fw=fx
        do iter=1,ITMAX
            xm=0.5_dp*(a+b)
            tol1=tol*abs(x)+ZEPS
            tol2=2.0_dp*tol1
            if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
                xmin=x
                brent=fx
                RETURN
            end if
            if (abs(e) > tol1) then
                r=(x-w)*(fx-fv)
                q=(x-v)*(fx-fw)
                p=(x-v)*q-(x-w)*r
                q=2.0_dp*(q-r)
                if (q > 0.0) p=-p
                q=abs(q)
                etemp=e
                e=d
                if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                    p <= q*(a-x) .or. p >= q*(b-x)) then
                    e=merge(a-x,b-x, x >= xm )
                    d=CGOLD*e
                else
                    d=p/q
                    u=x+d
                    if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
                end if
            else
                e=merge(a-x,b-x, x >= xm )
                d=CGOLD*e
            end if
            u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
            fu=func(u)
            if (fu <= fx) then
                if (u >= x) then
                    a=x
                else
                    b=x
                end if
                call shft(v,w,x,u)
                call shft(fv,fw,fx,fu)
            else
                if (u < x) then
                    a=u
                else
                    b=u
                end if
                if (fu <= fw .or. w == x) then
                    v=w
                    fv=fw
                    w=u
                    fw=fu
                else if (fu <= fv .or. v == x .or. v == w) then
                    v=u
                    fv=fu
                end if
            end if
        end do
        call nrerror('brent: exceed maximum iterations')
    CONTAINS
        !BL
        SUBROUTINE shft(a,b,c,d)
            REAL(dp), INTENT(OUT) :: a
            REAL(dp), INTENT(INOUT) :: b,c
            REAL(dp), INTENT(IN) :: d
            a=b
            b=c
            c=d
        END SUBROUTINE shft
    END FUNCTION brent

    SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
        USE nrtype; USE nrutil, ONLY : swap
        IMPLICIT NONE
        REAL(dp), INTENT(INOUT) :: ax,bx
        REAL(dp), INTENT(OUT) :: cx,fa,fb,fc
        PROCEDURE(template_function), POINTER, INTENT(IN) :: func
        REAL(dp), PARAMETER :: GOLD=1.618034_dp,GLIMIT=100.0_dp,TINY=1.0e-20_dp
        REAL(dp) :: fu,q,r,u,ulim
        fa=func(ax)
        fb=func(bx)
        if (fb > fa) then
            call swap(ax,bx)
            call swap(fa,fb)
        end if
        cx=bx+GOLD*(bx-ax)
        fc=func(cx)
        do
            if (fb < fc) RETURN
            r=(bx-ax)*(fb-fc)
            q=(bx-cx)*(fb-fa)
            u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_dp*sign(max(abs(q-r),TINY),q-r))
            ulim=bx+GLIMIT*(cx-bx)
            if ((bx-u)*(u-cx) > 0.0) then
                fu=func(u)
                if (fu < fc) then
                    ax=bx
                    fa=fb
                    bx=u
                    fb=fu
                    RETURN
                else if (fu > fb) then
                    cx=u
                    fc=fu
                    RETURN
                end if
                u=cx+GOLD*(cx-bx)
                fu=func(u)
            else if ((cx-u)*(u-ulim) > 0.0) then
                fu=func(u)
                if (fu < fc) then
                    bx=cx
                    cx=u
                    u=cx+GOLD*(cx-bx)
                    call shft(fb,fc,fu,func(u))
                end if
            else if ((u-ulim)*(ulim-cx) >= 0.0) then
                u=ulim
                fu=func(u)
            else
                u=cx+GOLD*(cx-bx)
                fu=func(u)
            end if
            call shft(ax,bx,cx,u)
            call shft(fa,fb,fc,fu)
        end do
    CONTAINS
        !BL
        SUBROUTINE shft(a,b,c,d)
            REAL(dp), INTENT(OUT) :: a
            REAL(dp), INTENT(INOUT) :: b,c
            REAL(dp), INTENT(IN) :: d
            a=b
            b=c
            c=d
        END SUBROUTINE shft
    END SUBROUTINE mnbrak

end module nr

MODULE  hmwk2
    use nrtype
    use nr
    implicit none
    INTEGER, PARAMETER :: N=74,LDA=N,LDB=N,NRHS=1
    REAL(DP), PARAMETER :: beta=1.0D0, beta2=0.96, eta=0.6, R=1.03D0, sigma=2.0D0

contains

    !-------------------------------------
    subroutine sub_mystop(calling)
        !-------------------------------------

        ! a personal stop subroutine. Makes it easier to edit behaviour of stop. All
        ! functions and subroutines should call this.
        !
        ! INPUTS: calling - a string indicating where this subroutine was called from

        CHARACTER (LEN=*), intent(in) :: calling
        print *, "STOP: ", calling
        STOP 0
    end subroutine sub_mystop

    !-------------------------------------
    SUBROUTINE q2aMethod1(func, w)
        !-------------------------------------

        PROCEDURE(template_function), POINTER, INTENT(IN) :: func
        REAL(DP), DIMENSION(N+2), INTENT(IN) :: w
        REAL(DP) :: fret
        INTEGER :: i
        REAL(DP) :: startPoint = 0.5D0, x
        REAL(DP), DIMENSION(N+2) :: assets, consumption

        print *, "-------------------------------------"
        print *, "forward"
        CALL mybrent(func, startPoint,fret)

        x=startPoint
        assets(1)=0
        assets(2)=x
        consumption(1)=R*assets(1)+w(1)-assets(2)

        open(unit=1,file="method1.csv")
        write(1,*) "age,assets,consumption,income"
        write(1,*) "16,0,",consumption(1),",",w(1)

        DO i=3,N+2
            assets(i)=assets(i-1)*(R+(beta*R)**(1/sigma))+w(i-1)-(beta*R)**(1/sigma)*(w(i-2)+R*assets(i-2))
            consumption(i-1)=R*assets(i-1)+w(i-1)-assets(i)
            write(1,*) 14+i,",",assets(i-1),",",consumption(i-1),",",w(i-1)
        END DO

        write(1,*) 14+i,",",assets(i-1),",",consumption(i-1),",",w(i-1)
        close(1)
    END SUBROUTINE q2aMethod1

    !-------------------------------------
    SUBROUTINE q2aMethod2(func, w)
        !-------------------------------------

        PROCEDURE(template_function), POINTER, INTENT(IN) :: func
        REAL(DP), DIMENSION(N+2), INTENT(IN) :: w
        REAL(DP) :: fret
        INTEGER :: i
        REAL(DP) :: startPoint = 0.5D0, x
        REAL(DP), DIMENSION(N+2) :: assets, consumption

        print *, "-------------------------------------"
        print *, "backward"
        CALL mybrent(func, startPoint,fret)

        x=startPoint
        assets(N+2)=0
        assets(N+1)=x

        open(unit=1,file="method2.csv")
        write(1,*) "age,assets,consumption,income"

        DO i=N+1,2,-1
            consumption(i)=R*assets(i)+w(i)-assets(i+1)
            assets(i-1)=1/R*(assets(i)-w(i-1)+(beta*R)**(-1/sigma)*(R*assets(i)+w(i)-assets(i+1)))
            write(1,*) 15+i,",",assets(i),",",consumption(i),",",w(i)
        END DO

        consumption(1)=R*assets(i)+w(1)-assets(2)
        write(1,*) "16,",assets(1),",",consumption(1),",",w(1)
        close(1)
    END SUBROUTINE q2aMethod2

    !-------------------------------------
    SUBROUTINE q2aMethod3(w)
        !-------------------------------------

        REAL(DP), DIMENSION(N+2), INTENT(IN) :: w
        REAL(DP), DIMENSION(LDA,N) :: A
        REAL(DP), DIMENSION(LDB,NRHS) :: b
        REAL(DP), DIMENSION(N) :: consumption
        INTEGER :: IPIV( N )
        INTEGER :: INFO,i,j

        print *, "-------------------------------------"
        print *, "Brute Force"

        DO i=1,N
            DO j=1,N
                if(i==j)then
                    A(j,i)=-1-(beta*R)**(-1/sigma)*R
                else if (j==(i+1))then
                    A(j,i)=R
                else if (j==(i-1)) then
                    A(j,i)=(beta*R)**(-1/sigma)
                else
                    A(j,i)=0.0D0
                end if
            END DO
            b(i,1)=(beta*R)**(-1/sigma)*w(i+1)-w(i)
        END DO

        CALL DGESV( N, NRHS, A, LDA, IPIV, b, LDB, INFO )

        IF( INFO.GT.0 ) THEN
            WRITE(*,*)'The diagonal element of the triangular factor of A,'
            WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
            WRITE(*,*)'A is singular; the solution could not be computed.'
            STOP
        END IF

        open(unit=1,file="method3.csv")
        write(1,*) "age,assets,consumption,income"

        consumption(1)=w(1)-b(1,1)
        DO i=1,N
            consumption(i+1)=R*b(i,1)+w(i+1)-b(i+1,1)
            write(1,*) 16+i,",",b(i,1),",",consumption(i+1),",",w(i+1)
        END DO
        close(1)
    END SUBROUTINE q2aMethod3

    !-------------------------------------
    SUBROUTINE q2b(w)
        !-------------------------------------

        REAL(DP), DIMENSION(N+2), INTENT(IN) :: w
        REAL(DP), DIMENSION(LDA,N) :: A
        REAL(DP), DIMENSION(LDB,NRHS) :: b
        REAL(DP), DIMENSION(N) :: Xv
        REAL(DP), DIMENSION(N+1) ::consumption,labour
        INTEGER :: IPIV( N )
        INTEGER :: INFO,i,j

        print *, "-------------------------------------"
        print *, "Part b"

        DO i=1,N
            Xv(i)=((beta2*R)*(w(i)/w(i+1))**((1-eta)*(1-sigma)))**(-1/sigma)
        END DO

        DO i=1,N
            DO j=1,N
                if(i==j)then
                    A(j,i)=-1-Xv(j)*R
                else if (j==(i+1))then
                    A(j,i)=R
                else if (j==(i-1)) then
                    A(j,i)=Xv(j)
                else
                    A(j,i)=0.0D0
                end if
            END DO
            b(i,1)=Xv(i)*w(i+1)-w(i)
        END DO

        CALL DGESV( N, NRHS, A, LDA, IPIV, b, LDB, INFO )

        IF( INFO.GT.0 ) THEN
            WRITE(*,*)'The diagonal element of the triangular factor of A,'
            WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
            WRITE(*,*)'A is singular; the solution could not be computed.'
            STOP
        END IF

        consumption(1)=eta*(w(1)-b(1,1))
        labour(1)=1-(1-eta)/eta*consumption(1)/w(1)

        DO i=1,N
            consumption(i+1)=eta*(R*b(i,1)+w(i+1)-b(i+1,1))
            labour(i+1)=1-(1-eta)/eta*consumption(i+1)/w(i+1)
        END DO

        open(unit=1,file="partb.csv")
        write(1,*) "age,assets,consumption,labour,income"
        write(1,*) "16,0,",consumption(1),",",labour(1),",",w(1)
        DO i=1,N
            write(1,*) 16+i,",",b(i,1),",",consumption(i+1),",",labour(i+1),",",w(i+1)
        END DO
        close(1)
    END SUBROUTINE q2b


    !-------------------------------------
    SUBROUTINE mybrent(func,startPoint,fret)
        !-------------------------------------

        !use brent
        !
        ! INPUTS:
        !   1) func - the function we are trying to find the min value of
        !   2) startPoint - the co-ordinates of the starting point
        PROCEDURE(template_function), POINTER, INTENT(in) :: func
        REAL(DP), INTENT(INOUT) :: startPoint
        REAL(DP), INTENT(OUT) :: fret
        REAL(DP):: tol=1.0e-12
        REAL(DP) :: temp1,temp2,temp3,f1,f2,f3

        temp1=startPoint
        temp2=0.5*startPoint
        call mnbrak(temp1, temp2, temp3, f1,f2,f3, func)
        fret=brent(temp1, temp2, temp3, func, tol, startPoint)

    END SUBROUTINE mybrent

END MODULE hmwk2

PROGRAM main
    use nrtype
    use hmwk2
    implicit none
    REAL(DP), DIMENSION(N+2) :: w
    PROCEDURE(template_function), POINTER :: func
    w =(/ 0.516_dp,0.560_dp,0.604_dp,0.648_dp,0.692_dp,0.736_dp,0.780_dp,&
        & 0.831_dp,0.883_dp,0.934_dp,0.986_dp,1.037_dp,1.089_dp,1.140_dp,&
        & 1.163_dp,1.186_dp,1.209_dp,1.232_dp,1.255_dp,1.278_dp,1.301_dp,&
        & 1.324_dp,1.347_dp,1.370_dp,1.372_dp,1.374_dp,1.376_dp,1.378_dp,&
        & 1.380_dp,1.382_dp,1.384_dp,1.386_dp,1.388_dp,1.390_dp,1.385_dp,&
        & 1.379_dp,1.374_dp,1.368_dp,1.363_dp,1.357_dp,1.352_dp,1.346_dp,&
        & 1.341_dp,1.335_dp,1.330_dp,1.304_dp,1.278_dp,1.252_dp,1.226_dp,&
        & 1.201_dp,1.175_dp,1.149_dp,1.123_dp,1.097_dp,1.071_dp,1.045_dp,&
        & 1.019_dp,0.994_dp,0.968_dp,0.942_dp,0.916_dp,0.890_dp,0.864_dp,&
        & 0.838_dp,0.812_dp,0.786_dp,0.761_dp,0.735_dp,0.709_dp,0.683_dp,&
        & 0.657_dp,0.631_dp,0.605_dp,0.579_dp,0.554_dp,0.000_dp /)

    func => forward
    CALL q2aMethod1(func, w)

    func => backward
    CALL q2aMethod2(func, w)

    CALL q2aMethod3(w)

    CALL q2b(w)

CONTAINS

    !-----------------
    FUNCTION forward(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: y - the value of the function at that point
        REAL(DP), INTENT(IN) :: point
        REAL(DP), DIMENSION(N+2) :: assets
        REAL(DP) :: x,z
        INTEGER :: i

        x=point
        assets(1)=0
        assets(2)=x

        DO i=3,N+2
            assets(i)=assets(i-1)*(R+(beta*R)**(1/sigma))+w(i-1)-(beta*R)**(1/sigma)*(w(i-2)+R*assets(i-2))
        END DO

        z=abs(assets(N+2))
    END FUNCTION forward

    !-----------------
    FUNCTION backward(point) RESULT(z)
        !
        ! the function we are trying to minimize
        ! INPUTS: point - the value at which we are evaluating the function
        ! OUTPUTS: y - the value of the function at that point
        REAL(DP), INTENT(IN) :: point
        REAL(DP), DIMENSION(N+2) :: assets
        REAL(DP) :: x,z
        INTEGER :: i

        x=point
        assets(N+2)=0
        assets(N+1)=x

        DO i=N,1,-1
            assets(i)=1/R*(assets(i+1)-w(i)+(beta*R)**(-1/sigma)*(R*assets(i+1)+w(i+1)-assets(i+2)))
        END DO

        z=abs(assets(1))
    END FUNCTION backward

end program main

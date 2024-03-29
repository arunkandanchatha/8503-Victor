MODULE nrutil
    USE nrtype
    IMPLICIT NONE
    INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
    INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
    INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
    INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
    INTEGER(I4B), PARAMETER :: NPAR_POLY=8
    INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8
    INTERFACE array_copy
        MODULE PROCEDURE array_copy_d, array_copy_i
    END INTERFACE
    INTERFACE swap
        MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
            swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
            masked_swap_rs,masked_swap_rv,masked_swap_rm
    END INTERFACE
    INTERFACE reallocate
        MODULE PROCEDURE reallocate_rv,reallocate_rm,&
            reallocate_iv,reallocate_im,reallocate_hv
    END INTERFACE
    INTERFACE imaxloc
        MODULE PROCEDURE imaxloc_r,imaxloc_i
    END INTERFACE
    INTERFACE assert
        MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
    END INTERFACE
    INTERFACE assert_eq
        MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
    END INTERFACE
    INTERFACE arth
        MODULE PROCEDURE arth_d, arth_i
    END INTERFACE
    INTERFACE geop
        MODULE PROCEDURE geop_d, geop_i, geop_c, geop_dv
    END INTERFACE
    INTERFACE cumsum
        MODULE PROCEDURE cumsum_r,cumsum_i
    END INTERFACE
    INTERFACE poly
        MODULE PROCEDURE poly_dd,poly_ddv,&
            poly_rc,poly_cc,poly_msk_ddv
    END INTERFACE
    INTERFACE poly_term
        MODULE PROCEDURE poly_term_rr,poly_term_cc
    END INTERFACE
    INTERFACE outerprod
        MODULE PROCEDURE outerprod_d
    END INTERFACE
    INTERFACE outerdiff
        MODULE PROCEDURE outerdiff_d,outerdiff_i
    END INTERFACE
    INTERFACE scatter_add
        MODULE PROCEDURE scatter_add_d
    END INTERFACE
    INTERFACE scatter_max
        MODULE PROCEDURE scatter_max_d
    END INTERFACE
    INTERFACE diagadd
        MODULE PROCEDURE diagadd_rv,diagadd_r
    END INTERFACE
    INTERFACE diagmult
        MODULE PROCEDURE diagmult_rv,diagmult_r
    END INTERFACE
    INTERFACE get_diag
        MODULE PROCEDURE get_diag_dv
    END INTERFACE
    INTERFACE put_diag
        MODULE PROCEDURE put_diag_rv, put_diag_r
    END INTERFACE
CONTAINS
    !BL
    SUBROUTINE array_copy_d(src,dest,n_copied,n_not_copied)
        REAL(DP), DIMENSION(:), INTENT(IN) :: src
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied=min(size(src),size(dest))
        n_not_copied=size(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
    END SUBROUTINE array_copy_d
    !BL
    SUBROUTINE array_copy_i(src,dest,n_copied,n_not_copied)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: src
        INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: dest
        INTEGER(I4B), INTENT(OUT) :: n_copied, n_not_copied
        n_copied=min(size(src),size(dest))
        n_not_copied=size(src)-n_copied
        dest(1:n_copied)=src(1:n_copied)
    END SUBROUTINE array_copy_i
    !BL
    !BL
    SUBROUTINE swap_i(a,b)
        INTEGER(I4B), INTENT(INOUT) :: a,b
        INTEGER(I4B) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_i
    !BL
    SUBROUTINE swap_r(a,b)
        REAL(DP), INTENT(INOUT) :: a,b
        REAL(DP) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_r
    !BL
    SUBROUTINE swap_rv(a,b)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
        REAL(DP), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_rv
    !BL
    SUBROUTINE swap_c(a,b)
        COMPLEX(SPC), INTENT(INOUT) :: a,b
        COMPLEX(SPC) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_c
    !BL
    SUBROUTINE swap_cv(a,b)
        COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_cv
    !BL
    SUBROUTINE swap_cm(a,b)
        COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_cm
    !BL
    SUBROUTINE swap_z(a,b)
        COMPLEX(DPC), INTENT(INOUT) :: a,b
        COMPLEX(DPC) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_z
    !BL
    SUBROUTINE swap_zv(a,b)
        COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_zv
    !BL
    SUBROUTINE swap_zm(a,b)
        COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
        COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
        dum=a
        a=b
        b=dum
    END SUBROUTINE swap_zm
    !BL
    SUBROUTINE masked_swap_rs(a,b,mask)
        REAL(DP), INTENT(INOUT) :: a,b
        LOGICAL(LGT), INTENT(IN) :: mask
        REAL(DP) :: swp
        if (mask) then
            swp=a
            a=b
            b=swp
        end if
    END SUBROUTINE masked_swap_rs
    !BL
    SUBROUTINE masked_swap_rv(a,b,mask)
        REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(DP), DIMENSION(size(a)) :: swp
        where (mask)
            swp=a
            a=b
            b=swp
        end where
    END SUBROUTINE masked_swap_rv
    !BL
    SUBROUTINE masked_swap_rm(a,b,mask)
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
        LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
        REAL(DP), DIMENSION(size(a,1),size(a,2)) :: swp
        where (mask)
            swp=a
            a=b
            b=swp
        end where
    END SUBROUTINE masked_swap_rm
    !BL
    !BL
    FUNCTION reallocate_rv(p,n)
        REAL(DP), DIMENSION(:), POINTER :: p, reallocate_rv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_rv(n),stat=ierr)
        if (ierr /= 0) call &
            nrerror('reallocate_rv: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_rv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    END FUNCTION reallocate_rv
    !BL
    FUNCTION reallocate_iv(p,n)
        INTEGER(I4B), DIMENSION(:), POINTER :: p, reallocate_iv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_iv(n),stat=ierr)
        if (ierr /= 0) call &
            nrerror('reallocate_iv: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    END FUNCTION reallocate_iv
    !BL
    FUNCTION reallocate_hv(p,n)
        CHARACTER(1), DIMENSION(:), POINTER :: p, reallocate_hv
        INTEGER(I4B), INTENT(IN) :: n
        INTEGER(I4B) :: nold,ierr
        allocate(reallocate_hv(n),stat=ierr)
        if (ierr /= 0) call &
            nrerror('reallocate_hv: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p)
        reallocate_hv(1:min(nold,n))=p(1:min(nold,n))
        deallocate(p)
    END FUNCTION reallocate_hv
    !BL
    FUNCTION reallocate_rm(p,n,m)
        REAL(DP), DIMENSION(:,:), POINTER :: p, reallocate_rm
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        allocate(reallocate_rm(n,m),stat=ierr)
        if (ierr /= 0) call &
            nrerror('reallocate_rm: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p,1)
        mold=size(p,2)
        reallocate_rm(1:min(nold,n),1:min(mold,m))=&
            p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    END FUNCTION reallocate_rm
    !BL
    FUNCTION reallocate_im(p,n,m)
        INTEGER(I4B), DIMENSION(:,:), POINTER :: p, reallocate_im
        INTEGER(I4B), INTENT(IN) :: n,m
        INTEGER(I4B) :: nold,mold,ierr
        allocate(reallocate_im(n,m),stat=ierr)
        if (ierr /= 0) call &
            nrerror('reallocate_im: problem in attempt to allocate memory')
        if (.not. associated(p)) RETURN
        nold=size(p,1)
        mold=size(p,2)
        reallocate_im(1:min(nold,n),1:min(mold,m))=&
            p(1:min(nold,n),1:min(mold,m))
        deallocate(p)
    END FUNCTION reallocate_im
    !BL
    FUNCTION ifirstloc(mask)
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        INTEGER(I4B) :: ifirstloc
        INTEGER(I4B), DIMENSION(1) :: loc
        loc=maxloc(merge(1,0,mask))
        ifirstloc=loc(1)
        if (.not. mask(ifirstloc)) ifirstloc=size(mask)+1
    END FUNCTION ifirstloc
    !BL
    FUNCTION imaxloc_r(arr)
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B) :: imaxloc_r
        INTEGER(I4B), DIMENSION(1) :: imax
        imax=maxloc(arr(:))
        imaxloc_r=imax(1)
    END FUNCTION imaxloc_r
    !BL
    FUNCTION imaxloc_i(iarr)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
        INTEGER(I4B), DIMENSION(1) :: imax
        INTEGER(I4B) :: imaxloc_i
        imax=maxloc(iarr(:))
        imaxloc_i=imax(1)
    END FUNCTION imaxloc_i
    !BL
    FUNCTION iminloc(arr)
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), DIMENSION(1) :: imin
        INTEGER(I4B) :: iminloc
        imin=minloc(arr(:))
        iminloc=imin(1)
    END FUNCTION iminloc
    !BL
    SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        if (.not. n1) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                string
            STOP 'program terminated by assert1'
        end if
    END SUBROUTINE assert1
    !BL
    SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        if (.not. (n1 .and. n2)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                string
            STOP 'program terminated by assert2'
        end if
    END SUBROUTINE assert2
    !BL
    SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        if (.not. (n1 .and. n2 .and. n3)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                string
            STOP 'program terminated by assert3'
        end if
    END SUBROUTINE assert3
    !BL
    SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                string
            STOP 'program terminated by assert4'
        end if
    END SUBROUTINE assert4
    !BL
    SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        if (.not. all(n)) then
            write (*,*) 'nrerror: an assertion failed with this tag:', &
                string
            STOP 'program terminated by assert_v'
        end if
    END SUBROUTINE assert_v
    !BL
    FUNCTION assert_eq2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2
        INTEGER :: assert_eq2
        if (n1 == n2) then
            assert_eq2=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                string
            STOP 'program terminated by assert_eq2'
        end if
    END FUNCTION assert_eq2
    !BL
    FUNCTION assert_eq3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3
        INTEGER :: assert_eq3
        if (n1 == n2 .and. n2 == n3) then
            assert_eq3=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                string
            STOP 'program terminated by assert_eq3'
        end if
    END FUNCTION assert_eq3
    !BL
    FUNCTION assert_eq4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, INTENT(IN) :: n1,n2,n3,n4
        INTEGER :: assert_eq4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq4=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                string
            STOP 'program terminated by assert_eq4'
        end if
    END FUNCTION assert_eq4
    !BL
    FUNCTION assert_eqn(nn,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        INTEGER, DIMENSION(:), INTENT(IN) :: nn
        INTEGER :: assert_eqn
        if (all(nn(2:) == nn(1))) then
            assert_eqn=nn(1)
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag:', &
                string
            STOP 'program terminated by assert_eqn'
        end if
    END FUNCTION assert_eqn
    !BL
    SUBROUTINE nrerror(string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        write (*,*) 'nrerror: ',string
        print*, 'program terminated by nrerror'
        stop 0
    END SUBROUTINE nrerror
    !BL
    FUNCTION arth_d(first,increment,n)
        REAL(DP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: arth_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        if (n > 0) arth_d(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
                arth_d(k)=arth_d(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
                arth_d(k)=arth_d(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                if (k >= n) exit
                k2=k+k
                arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
                temp=temp+temp
                k=k2
            end do
        end if
    END FUNCTION arth_d
    !BL
    FUNCTION arth_i(first,increment,n)
        INTEGER(I4B), INTENT(IN) :: first,increment,n
        INTEGER(I4B), DIMENSION(n) :: arth_i
        INTEGER(I4B) :: k,k2,temp
        if (n > 0) arth_i(1)=first
        if (n <= NPAR_ARTH) then
            do k=2,n
                arth_i(k)=arth_i(k-1)+increment
            end do
        else
            do k=2,NPAR2_ARTH
                arth_i(k)=arth_i(k-1)+increment
            end do
            temp=increment*NPAR2_ARTH
            k=NPAR2_ARTH
            do
                if (k >= n) exit
                k2=k+k
                arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
                temp=temp+temp
                k=k2
            end do
        end if
    END FUNCTION arth_i
    !BL
    FUNCTION geop_d(first,factor,n)
#if 0
        REAL(qp), INTENT(IN) :: first,factor
#else
        REAL(DP), INTENT(IN) :: first,factor
#endif
        INTEGER(I4B), INTENT(IN) :: n
#if 0
        REAL(qp), DIMENSION(n) :: geop_d
#else
        REAL(DP), DIMENSION(n) :: geop_d
#endif
        INTEGER(I4B) :: k,k2
#if 0
        REAL(qp) :: temp
#else
        REAL(DP) :: temp
#endif
        if (n > 0) geop_d(1)=first
        if (n <= NPAR_GEOP) then
            do k=2,n
                geop_d(k)=geop_d(k-1)*factor
            end do
        else
            do k=2,NPAR2_GEOP
                geop_d(k)=geop_d(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                if (k >= n) exit
                k2=k+k
                geop_d(k+1:min(k2,n))=temp*geop_d(1:min(k,n-k))
                temp=temp*temp
                k=k2
            end do
        end if
    END FUNCTION geop_d
    !BL
    FUNCTION geop_i(first,factor,n)
        INTEGER(I4B), INTENT(IN) :: first,factor,n
        INTEGER(I4B), DIMENSION(n) :: geop_i
        INTEGER(I4B) :: k,k2,temp
        if (n > 0) geop_i(1)=first
        if (n <= NPAR_GEOP) then
            do k=2,n
                geop_i(k)=geop_i(k-1)*factor
            end do
        else
            do k=2,NPAR2_GEOP
                geop_i(k)=geop_i(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                if (k >= n) exit
                k2=k+k
                geop_i(k+1:min(k2,n))=temp*geop_i(1:min(k,n-k))
                temp=temp*temp
                k=k2
            end do
        end if
    END FUNCTION geop_i
    !BL
    FUNCTION geop_c(first,factor,n)
#if 0
        COMPLEX(dp), INTENT(IN) :: first,factor
#else
        COMPLEX(DP), INTENT(IN) :: first,factor
#endif
        INTEGER(I4B), INTENT(IN) :: n
#if 0
        COMPLEX(dp), DIMENSION(n) :: geop_c
#else
        COMPLEX(DP), DIMENSION(n) :: geop_c
#endif
        INTEGER(I4B) :: k,k2
#if 0
        COMPLEX(dp) :: temp
#else
        COMPLEX(DP) :: temp
#endif
        if (n > 0) geop_c(1)=first
        if (n <= NPAR_GEOP) then
            do k=2,n
                geop_c(k)=geop_c(k-1)*factor
            end do
        else
            do k=2,NPAR2_GEOP
                geop_c(k)=geop_c(k-1)*factor
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                if (k >= n) exit
                k2=k+k
                geop_c(k+1:min(k2,n))=temp*geop_c(1:min(k,n-k))
                temp=temp*temp
                k=k2
            end do
        end if
    END FUNCTION geop_c
    !BL
    FUNCTION geop_dv(first,factor,n)
#if 0
        REAL(qp), DIMENSION(:), INTENT(IN) :: first,factor
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: first,factor
#endif
        INTEGER(I4B), INTENT(IN) :: n
#if 0
        REAL(qp), DIMENSION(size(first),n) :: geop_dv
#else
        REAL(DP), DIMENSION(size(first),n) :: geop_dv
#endif
        INTEGER(I4B) :: k,k2
#if 0
        REAL(qp), DIMENSION(size(first)) :: temp
#else
        REAL(DP), DIMENSION(size(first)) :: temp
#endif
        if (n > 0) geop_dv(:,1)=first(:)
        if (n <= NPAR_GEOP) then
            do k=2,n
                geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
            end do
        else
            do k=2,NPAR2_GEOP
                geop_dv(:,k)=geop_dv(:,k-1)*factor(:)
            end do
            temp=factor**NPAR2_GEOP
            k=NPAR2_GEOP
            do
                if (k >= n) exit
                k2=k+k
                geop_dv(:,k+1:min(k2,n))=geop_dv(:,1:min(k,n-k))*&
                    spread(temp,2,size(geop_dv(:,1:min(k,n-k)),2))
                temp=temp*temp
                k=k2
            end do
        end if
    END FUNCTION geop_dv
    !BL
    !BL
    RECURSIVE FUNCTION cumsum_r(arr,seed) RESULT(ans)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: arr
        REAL(dp), OPTIONAL, INTENT(IN) :: seed
        REAL(dp), DIMENSION(size(arr)) :: ans
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        REAL(DP), OPTIONAL, INTENT(IN) :: seed
        REAL(DP), DIMENSION(size(arr)) :: ans
#endif
        INTEGER(I4B) :: n,j
#if 0
        REAL(dp) :: sd
#else
        REAL(DP) :: sd
#endif
        n=size(arr)
        if (n == 0_i4b) RETURN
#if 0
        sd=0.0_dp
#else
        sd=0.0_dp
#endif
        if (present(seed)) sd=seed
        ans(1)=arr(1)+sd
        if (n < NPAR_CUMSUM) then
            do j=2,n
                ans(j)=ans(j-1)+arr(j)
            end do
        else
            ans(2:n:2)=cumsum_r(arr(2:n:2)+arr(1:n-1:2),sd)
            ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
        end if
    END FUNCTION cumsum_r
    !BL
    RECURSIVE FUNCTION cumsum_i(arr,seed) RESULT(ans)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
        INTEGER(I4B), DIMENSION(size(arr)) :: ans
        INTEGER(I4B) :: n,j,sd
        n=size(arr)
        if (n == 0_i4b) RETURN
        sd=0_i4b
        if (present(seed)) sd=seed
        ans(1)=arr(1)+sd
        if (n < NPAR_CUMSUM) then
            do j=2,n
                ans(j)=ans(j-1)+arr(j)
            end do
        else
            ans(2:n:2)=cumsum_i(arr(2:n:2)+arr(1:n-1:2),sd)
            ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
        end if
    END FUNCTION cumsum_i
    !BL
    !BL
    RECURSIVE FUNCTION cumprod(arr,seed) RESULT(ans)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: arr
        REAL(dp), OPTIONAL, INTENT(IN) :: seed
        REAL(dp), DIMENSION(size(arr)) :: ans
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: arr
        REAL(DP), OPTIONAL, INTENT(IN) :: seed
        REAL(DP), DIMENSION(size(arr)) :: ans
#endif
        INTEGER(I4B) :: n,j
#if 0
        REAL(dp) :: sd
#else
        REAL(DP) :: sd
#endif
        n=size(arr)
        if (n == 0_i4b) RETURN
#if 0
        sd=1.0_dp
#else
        sd=1.0_dp
#endif
        if (present(seed)) sd=seed
        ans(1)=arr(1)*sd
        if (n < NPAR_CUMPROD) then
            do j=2,n
                ans(j)=ans(j-1)*arr(j)
            end do
        else
            ans(2:n:2)=cumprod(arr(2:n:2)*arr(1:n-1:2),sd)
            ans(3:n:2)=ans(2:n-1:2)*arr(3:n:2)
        end if
    END FUNCTION cumprod
    !BL
    FUNCTION poly_dd(x,coeffs)
#if 0
        REAL(qp), INTENT(IN) :: x
        REAL(qp), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(qp) :: poly_dd
        REAL(qp) :: pow
        REAL(qp), DIMENSION(:), ALLOCATABLE :: vec
#else
        REAL(DP), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(DP) :: poly_dd
        REAL(DP) :: pow
        REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
#endif
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
#if 0
            poly_dd=0.0_qp
#else
            poly_dd=0.0_dp
#endif
        else if (n < NPAR_POLY) then
            poly_dd=coeffs(n)
            do i=n-1,1,-1
                poly_dd=x*poly_dd+coeffs(i)
            end do
        else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
#if 0
                vec(n+1)=0.0_qp
#else
                vec(n+1)=0.0_dp
#endif
                nn=ishft(n+1,-1)
                vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                if (nn == 1) exit
                pow=pow*pow
                n=nn
            end do
            poly_dd=vec(1)
            deallocate(vec)
        end if
    END FUNCTION poly_dd
    !BL
    FUNCTION poly_rc(x,coeffs)
#if 0
        COMPLEX(dpc), INTENT(IN) :: x
        REAL(dp), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(dpc) :: poly_rc
        COMPLEX(dpc) :: pow
        COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: vec
#else
        COMPLEX(SPC), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_rc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
#endif
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
#if 0
            poly_rc=0.0_dp
#else
            poly_rc=0.0_dp
#endif
        else if (n < NPAR_POLY) then
            poly_rc=coeffs(n)
            do i=n-1,1,-1
                poly_rc=x*poly_rc+coeffs(i)
            end do
        else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
#if 0
                vec(n+1)=0.0_dp
#else
                vec(n+1)=0.0_dp
#endif
                nn=ishft(n+1,-1)
                vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                if (nn == 1) exit
                pow=pow*pow
                n=nn
            end do
            poly_rc=vec(1)
            deallocate(vec)
        end if
    END FUNCTION poly_rc
    !BL
    FUNCTION poly_cc(x,coeffs)
#if 0
        COMPLEX(dpc), INTENT(IN) :: x
        COMPLEX(dpc), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(dpc) :: poly_cc
        COMPLEX(dpc) :: pow
        COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: vec
#else
        COMPLEX(SPC), INTENT(IN) :: x
        COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_cc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
#endif
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
#if 0
            poly_cc=0.0_dp
#else
            poly_cc=0.0_dp
#endif
        else if (n < NPAR_POLY) then
            poly_cc=coeffs(n)
            do i=n-1,1,-1
                poly_cc=x*poly_cc+coeffs(i)
            end do
        else
            allocate(vec(n+1))
            pow=x
            vec(1:n)=coeffs
            do
#if 0
                vec(n+1)=0.0_dp
#else
                vec(n+1)=0.0_dp
#endif
                nn=ishft(n+1,-1)
                vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                if (nn == 1) exit
                pow=pow*pow
                n=nn
            end do
            poly_cc=vec(1)
            deallocate(vec)
        end if
    END FUNCTION poly_cc
    !BL
    FUNCTION poly_ddv(x,coeffs)
#if 0
        REAL(qp), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(qp), DIMENSION(size(x)) :: poly_ddv
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(DP), DIMENSION(size(x)) :: poly_ddv
#endif
        INTEGER(I4B) :: i,n,m
        m=size(coeffs)
        n=size(x)
        if (m <= 0) then
#if 0
            poly_ddv=0.0_qp
#else
            poly_ddv=0.0_dp
#endif
        else if (m < n .or. m < NPAR_POLY) then
            poly_ddv=coeffs(m)
            do i=m-1,1,-1
                poly_ddv=x*poly_ddv+coeffs(i)
            end do
        else
            do i=1,n
                poly_ddv(i)=poly_dd(x(i),coeffs)
            end do
        end if
    END FUNCTION poly_ddv
    !BL
    FUNCTION poly_msk_ddv(x,coeffs,mask)
#if 0
        REAL(qp), DIMENSION(:), INTENT(IN) :: coeffs,x
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
#endif
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
#if 0
        REAL(qp), DIMENSION(size(x)) :: poly_msk_ddv
        poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_qp)
#else
        REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
        poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
#endif
    END FUNCTION poly_msk_ddv
    !BL
    !BL
    RECURSIVE FUNCTION poly_term_rr(a,b) RESULT(u)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: a
        REAL(dp), INTENT(IN) :: b
        REAL(dp), DIMENSION(size(a)) :: u
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: a
        REAL(DP), INTENT(IN) :: b
        REAL(DP), DIMENSION(size(a)) :: u
#endif
        INTEGER(I4B) :: n,j
        n=size(a)
        if (n <= 0) RETURN
        u(1)=a(1)
        if (n < NPAR_POLYTERM) then
            do j=2,n
                u(j)=a(j)+b*u(j-1)
            end do
        else
            u(2:n:2)=poly_term_rr(a(2:n:2)+a(1:n-1:2)*b,b*b)
            u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
        end if
    END FUNCTION poly_term_rr
    !BL
    RECURSIVE FUNCTION poly_term_cc(a,b) RESULT(u)
#if 0
        COMPLEX(dpc), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(dpc), INTENT(IN) :: b
        COMPLEX(dpc), DIMENSION(size(a)) :: u
#else
        COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: a
        COMPLEX(SPC), INTENT(IN) :: b
        COMPLEX(SPC), DIMENSION(size(a)) :: u
#endif
        INTEGER(I4B) :: n,j
        n=size(a)
        if (n <= 0) RETURN
        u(1)=a(1)
        if (n < NPAR_POLYTERM) then
            do j=2,n
                u(j)=a(j)+b*u(j-1)
            end do
        else
            u(2:n:2)=poly_term_cc(a(2:n:2)+a(1:n-1:2)*b,b*b)
            u(3:n:2)=a(3:n:2)+b*u(2:n-1:2)
        end if
    END FUNCTION poly_term_cc
    !BL
    !BL
    FUNCTION zroots_unity(n,nn)
        INTEGER(I4B), INTENT(IN) :: n,nn
#if 0
        COMPLEX(dpc), DIMENSION(nn) :: zroots_unity
#else
        COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
#endif
        INTEGER(I4B) :: k
#if 0
        REAL(dp) :: theta
#else
        REAL(DP) :: theta
#endif
        zroots_unity(1)=1.0
        theta=TWOPI/n
        k=1
        do
            if (k >= nn) exit
#if 0
            zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),dpc)
#else
            zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
#endif
            zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
                zroots_unity(2:min(k,nn-k))
            k=2*k
        end do
    END FUNCTION zroots_unity
    !BL
    FUNCTION outerprod_d(a,b)
#if 0
        REAL(qp), DIMENSION(:), INTENT(IN) :: a,b
        REAL(qp), DIMENSION(size(a),size(b)) :: outerprod_d
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
#endif
        outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
            spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerprod_d
    !BL
    FUNCTION outerdiv(a,b)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: a,b
        REAL(dp), DIMENSION(size(a),size(b)) :: outerdiv
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(size(a),size(b)) :: outerdiv
#endif
        outerdiv = spread(a,dim=2,ncopies=size(b)) / &
            spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerdiv
    !BL
    FUNCTION outersum(a,b)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: a,b
        REAL(dp), DIMENSION(size(a),size(b)) :: outersum
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(size(a),size(b)) :: outersum
#endif
        outersum = spread(a,dim=2,ncopies=size(b)) + &
            spread(b,dim=1,ncopies=size(a))
    END FUNCTION outersum
    !BL
    FUNCTION outerdiff_d(a,b)
#if 0
        REAL(qp), DIMENSION(:), INTENT(IN) :: a,b
        REAL(qp), DIMENSION(size(a),size(b)) :: outerdiff_d
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
        REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
#endif
        outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
            spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerdiff_d
    !BL
    FUNCTION outerdiff_i(a,b)
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
        INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
        outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
            spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerdiff_i
    !BL
    FUNCTION outerand(a,b)
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
        LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
        outerand = spread(a,dim=2,ncopies=size(b)) .and. &
            spread(b,dim=1,ncopies=size(a))
    END FUNCTION outerand
    !BL
    SUBROUTINE scatter_add_d(dest,source,dest_index)
#if 0
        REAL(qp), DIMENSION(:), INTENT(OUT) :: dest
        REAL(qp), DIMENSION(:), INTENT(IN) :: source
#else
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(DP), DIMENSION(:), INTENT(IN) :: source
#endif
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n=assert_eq2(size(source),size(dest_index),'scatter_add_d')
        m=size(dest)
        do j=1,n
            i=dest_index(j)
            if (i > 0 .and. i <= m) dest(i)=dest(i)+source(j)
        end do
    END SUBROUTINE scatter_add_d
    SUBROUTINE scatter_max_d(dest,source,dest_index)
#if 0
        REAL(qp), DIMENSION(:), INTENT(OUT) :: dest
        REAL(qp), DIMENSION(:), INTENT(IN) :: source
#else
        REAL(DP), DIMENSION(:), INTENT(OUT) :: dest
        REAL(DP), DIMENSION(:), INTENT(IN) :: source
#endif
        INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dest_index
        INTEGER(I4B) :: m,n,j,i
        n=assert_eq2(size(source),size(dest_index),'scatter_max_d')
        m=size(dest)
        do j=1,n
            i=dest_index(j)
            if (i > 0 .and. i <= m) dest(i)=max(dest(i),source(j))
        end do
    END SUBROUTINE scatter_max_d
    !BL
    SUBROUTINE diagadd_rv(mat,diag)
#if 0
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(dp), DIMENSION(:), INTENT(IN) :: diag
#else
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(DP), DIMENSION(:), INTENT(IN) :: diag
#endif
        INTEGER(I4B) :: j,n
        n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagadd_rv')
        do j=1,n
            mat(j,j)=mat(j,j)+diag(j)
        end do
    END SUBROUTINE diagadd_rv
    !BL
    SUBROUTINE diagadd_r(mat,diag)
#if 0
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(dp), INTENT(IN) :: diag
#else
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(DP), INTENT(IN) :: diag
#endif
        INTEGER(I4B) :: j,n
        n = min(size(mat,1),size(mat,2))
        do j=1,n
            mat(j,j)=mat(j,j)+diag
        end do
    END SUBROUTINE diagadd_r
    !BL
    SUBROUTINE diagmult_rv(mat,diag)
#if 0
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(dp), DIMENSION(:), INTENT(IN) :: diag
#else
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(DP), DIMENSION(:), INTENT(IN) :: diag
#endif
        INTEGER(I4B) :: j,n
        n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
        do j=1,n
            mat(j,j)=mat(j,j)*diag(j)
        end do
    END SUBROUTINE diagmult_rv
    !BL
    SUBROUTINE diagmult_r(mat,diag)
#if 0
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(dp), INTENT(IN) :: diag
#else
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
        REAL(DP), INTENT(IN) :: diag
#endif
        INTEGER(I4B) :: j,n
        n = min(size(mat,1),size(mat,2))
        do j=1,n
            mat(j,j)=mat(j,j)*diag
        end do
    END SUBROUTINE diagmult_r
    !BL
    FUNCTION get_diag_dv(mat)
#if 0
        REAL(qp), DIMENSION(:,:), INTENT(IN) :: mat
        REAL(qp), DIMENSION(size(mat,1)) :: get_diag_dv
#else
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
        REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
#endif
        INTEGER(I4B) :: j
        j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
        do j=1,size(mat,1)
            get_diag_dv(j)=mat(j,j)
        end do
    END FUNCTION get_diag_dv
    !BL
    SUBROUTINE put_diag_rv(diagv,mat)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: diagv
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: mat
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: diagv
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
#endif
        INTEGER(I4B) :: j,n
        n=assert_eq2(size(diagv),min(size(mat,1),size(mat,2)),'put_diag_rv')
        do j=1,n
            mat(j,j)=diagv(j)
        end do
    END SUBROUTINE put_diag_rv
    !BL
    SUBROUTINE put_diag_r(scal,mat)
#if 0
        REAL(dp), INTENT(IN) :: scal
        REAL(dp), DIMENSION(:,:), INTENT(INOUT) :: mat
#else
        REAL(DP), INTENT(IN) :: scal
        REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: mat
#endif
        INTEGER(I4B) :: j,n
        n = min(size(mat,1),size(mat,2))
        do j=1,n
            mat(j,j)=scal
        end do
    END SUBROUTINE put_diag_r
    !BL
    SUBROUTINE unit_matrix(mat)
#if 0
        REAL(dp), DIMENSION(:,:), INTENT(OUT) :: mat
#else
        REAL(DP), DIMENSION(:,:), INTENT(OUT) :: mat
#endif
        INTEGER(I4B) :: i,n
        n=min(size(mat,1),size(mat,2))
#if 0
        mat(:,:)=0.0_dp
#else
        mat(:,:)=0.0_dp
#endif
        do i=1,n
#if 0
            mat(i,i)=1.0_dp
#else
            mat(i,i)=1.0_dp
#endif
        end do
    END SUBROUTINE unit_matrix
    !BL
    FUNCTION upper_triangle(j,k,extra)
        INTEGER(I4B), INTENT(IN) :: j,k
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
        LOGICAL(LGT), DIMENSION(j,k) :: upper_triangle
        INTEGER(I4B) :: n
        n=0
        if (present(extra)) n=extra
        upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
    END FUNCTION upper_triangle
    !BL
    FUNCTION lower_triangle(j,k,extra)
        INTEGER(I4B), INTENT(IN) :: j,k
        INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
        LOGICAL(LGT), DIMENSION(j,k) :: lower_triangle
        INTEGER(I4B) :: n
        n=0
        if (present(extra)) n=extra
        lower_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) > -n)
    END FUNCTION lower_triangle
    !BL
    FUNCTION vabs(v)
#if 0
        REAL(dp), DIMENSION(:), INTENT(IN) :: v
        REAL(dp) :: vabs
#else
        REAL(DP), DIMENSION(:), INTENT(IN) :: v
        REAL(DP) :: vabs
#endif
        vabs=sqrt(dot_product(v,v))
    END FUNCTION vabs
!BL
END MODULE nrutil

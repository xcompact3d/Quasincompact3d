    subroutine test_min_max(name,text,array_tmp,i_size_array_tmp)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: ierror, i, i_size_array_tmp
    real(mytype) :: max_tmp, min_tmp, tot_tmp
    real(mytype), dimension(i_size_array_tmp) :: array_tmp
    character(len=5) :: name
    character(len=15) :: text

    max_tmp=-0.000000000000000001_mytype
    tot_tmp=0._mytype
    min_tmp=+1000000000000000000._mytype

    write(*,*) name,size(array_tmp)

    do i=1,size(array_tmp)
       max_tmp=max(max_tmp,array_tmp(i))
       tot_tmp=tot_tmp + array_tmp(i)
       min_tmp=min(min_tmp,array_tmp(i))
    enddo
    write(*,*) trim(text)//' Max Tot Min ',name,max_tmp,tot_tmp,min_tmp,nrank
    call flush(6)

    return
    end subroutine test_min_max

    subroutine test_min_max_comp(name,text,array_tmp,i_size_array_tmp)

    USE param
    USE variables
    USE decomp_2d
    USE MPI

    implicit none

    integer :: ierror, i, i_size_array_tmp
    real(mytype) :: tmp_max_real, tmp_max_complex, tmp_min_real, tmp_min_complex, tot_tmp_real, tot_tmp_complex
    complex(mytype) :: max_tmp, tot_tmp, min_tmp 
    complex(mytype), dimension(i_size_array_tmp) :: array_tmp
    character(len=5) :: name
    character(len=15) :: text

    tmp_max_real=-0.000000000000000001_mytype
    tmp_max_complex=-0.000000000000000001_mytype
    max_tmp=cmplx(tmp_max_real, tmp_max_complex)
    tot_tmp_real=0._mytype
    tot_tmp_complex=0._mytype
    tot_tmp=cmplx(tot_tmp_real, tot_tmp_complex)
    tmp_min_real=+1000000000000000000._mytype
    tmp_min_complex=+1000000000000000000._mytype
    min_tmp=cmplx(tmp_min_real, tmp_min_complex)

    do i=1,size(array_tmp)
       tmp_max_real=max(tmp_max_real, real(array_tmp(i)))
       tmp_max_complex=max(tmp_max_complex, aimag(array_tmp(i)))
       tot_tmp_real=tot_tmp_real + real(array_tmp(i))
       tot_tmp_complex=tot_tmp_complex + aimag(array_tmp(i))
       tmp_min_real=min(tmp_min_real, real(array_tmp(i)))
       tmp_min_complex=min(tmp_min_complex, aimag(array_tmp(i)))
    enddo
    max_tmp=cmplx(tmp_max_real, tmp_max_complex)
    tot_tmp=cmplx(tot_tmp_real, tot_tmp_complex)
    min_tmp=cmplx(tmp_min_real, tmp_min_complex)
    write(*,*) trim(text)//' Max Tot Min ',name,max_tmp,tot_tmp,min_tmp,nrank
    call flush(6)

    return
    end subroutine test_min_max_comp

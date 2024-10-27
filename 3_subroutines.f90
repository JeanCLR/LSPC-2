!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++       S U B R O U T I N E S        ++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++    Calculo de los Coeficientes de Fourier    +++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Parte que implementa el método numérico.
! mp : m_max + 1
!
subroutine get_coef(n,mp,k,coef)
    use constants
    implicit none
    integer(dp) :: i,j
    integer(dp),intent(in) :: n,mp
    real(dp) :: error,w
    real(dp),dimension(n),intent(inout) :: k
    complex(dp),dimension(mp,2*n),intent(inout) :: coef
    real(dp),dimension(:,:),allocatable :: A
    real(dp),dimension(:),allocatable :: sol_aux1,B
    
    error = (10.0_dp**(-6.0_dp))
    w = 1.5_dp

    allocate(A(n*4,n*4))
    allocate(B(n*4))
    allocate(sol_aux1(n*4))

    
    call m_zeros(A,n*4,n*4)
    call v_zeros(B,n*4)
    call m_zeros_complex(coef,mp,2*n)
    
    do i=0,mp-1,1
        call matrix_problem(A,B,k,n,i)
        call sor(A,B,n*4,error,w,.true.,sol_aux1)
        do j=1,2*n,1
            coef(i+1,j)=complex(sol_aux1(j),sol_aux1(j+(2*n)))
        end do
    end do

    deallocate(A)
    deallocate(B)
    deallocate(sol_aux1)

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++    Vector de puntos igualmente espaciados    +++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Se incluyen los extremos 'from' y 'to'.
!
subroutine linspace(from,to,vec,n)
    use constants
    implicit none
    integer(dp) :: i
    integer(dp), intent(in) :: n
    real(dp), intent(in) :: from, to
    real(dp), dimension(n), intent(inout) :: vec
    real(dp) :: step

    select case (n)
        case (0)
            return
        case (1)
            vec(1) = from
            return
    end select

    step = (to - from) / real(n-1, dp)

    do i=1,n,1
        vec(i) = from + (step*real(i-1, dp))
    end do

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++    Completa el Sistema a Resolver    +++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Depende del modo (m) que corresponde resolver
!
subroutine matrix_problem(A,B,k,n,m)
    use constants
    implicit none
    integer(dp) :: i, j
    integer(dp),intent(in) :: n, m
    real(dp),dimension(0:(n*4)-1,0:(n*4)-1),intent(inout) :: A
    real(dp),dimension(0:(n*4)-1),intent(inout) :: B
    real(dp),dimension(n),intent(inout) :: k
    real(dp) :: alpha, i_re, aux1, aux2, aux3
    complex(dp) :: aux4, contorno_simpsom
    
    alpha=0.0_dp
    if (mod(m,4_dp)==0_dp) then
        alpha=4.0_dp
    end if
    
    A(0,0) = 4.0_dp
    A(2*n,2*n) = 4.0_dp
    A(0,1) = (-1.0_dp)*alpha
    A(2*n,(2*n)+1) = (-1.0_dp)*alpha

    do i=1,(2*n)-1,1
        if (mod(i,2)/=0) then
            i_re = real(i, dp)
            aux1 = ((2.0_dp*i_re)-1.0_dp)/2.0_dp
            aux2 = ((-2.0_dp*(i_re**2.0_dp))-(real(m, dp)**2.0_dp))/i_re
            aux3 = ((2.0_dp*i_re)+1.0_dp)/2.0_dp
            A(i,i-1) = aux1
            A(i+(2*n),i+(2*n)-1)= aux1
            A(i,i) = aux2
            A(i+(2*n),i+(2*n))= aux2
            if (i/=(2*n)-1) then
                A(i,i+1) = aux3
                A(i+(2*n),i+(2*n)+1)= aux3
            end if
            if (i==(2*n)-1) then
                aux4 = contorno_simpsom(m)
                B(i) = (-1.0_dp)*aux3*real(aux4, dp)
                B(i+(2*n)) = (-1.0_dp)*aux3*aimag(aux4)
            end if
        end if
        if (mod(i,2)==0) then
            j = int(i/2)
            A(i,i-2) = k(j)
            A(i+(2*n),i+(2*n)-2)= k(j)
            A(i,i-1) = (-4.0_dp)*k(j)
            A(i+(2*n),i+(2*n)-1)= (-4.0_dp)*k(j)
            A(i,i) = 3.0_dp*(k(j)+k(j+1))
            A(i+(2*n),i+(2*n))= 3.0_dp*(k(j)+k(j+1))
            A(i,i+1) = (-4.0_dp)*k(j+1)
            A(i+(2*n),i+(2*n)+1)= (-4.0_dp)*k(j+1)
            if (i/=(2*n)-2) then
                A(i,i+2) = k(j+1)
                A(i+(2*n),i+(2*n)+2)= k(j+1)
            end if
            if (i==(2*n)-2) then
                aux4 = contorno_simpsom(m)
                B(i) = (-1.0_dp)*k(j+1)*real(aux4, dp)
                B(i+(2*n)) = (-1.0_dp)*k(j+1)*aimag(aux4)
            end if
        end if
    end do

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++        M e t o d o    S O R        ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
subroutine sor(A,B,n,error,w,logic,sol)
    use constants
    implicit none
    integer(dp),intent(in) :: n
    integer(dp) :: i,j,ite
    real(dp),intent(in) :: error,  w
    real(dp) :: delta,cuad,suma1,suma2,suma3
    real(dp),dimension(n,n),intent(in) :: A
    real(dp),dimension(n),intent(in) :: B
    real(dp),dimension(n),intent(inout) :: sol
    real(dp),dimension(n) :: aux,aux1,dif2
    logical,intent(in) :: logic
    ite=0 ; delta=100.0_dp

    call v_zeros(sol,n)
    call v_zeros(aux,n)
    call v_zeros(aux1,n)
    call v_zeros(dif2,n)

    do while (delta>=error)
        aux1=sol
        do i=1,n,1
            suma1 = 0.0_dp
            do j=1,i-1,1
                suma1 = suma1 + (A(i,j)*sol(j))
            end do
            suma2 = 0.0_dp
            do j=i+1,n,1
                suma2 = suma2 + (A(i,j)*sol(j))
            end do
            suma3 = 0.0_dp
            do j=1,i-1,1
                suma3 = suma3 + (A(i,j)*aux(j))
            end do
            aux(i)=(B(i)-((1.0_dp-w)*suma1)-suma2-(w*suma3))/A(i,i)
        end do
        sol=aux
        cuad = 0.0_dp
        dif2=(sol-aux1)
        do i=1,n,1
            cuad = cuad + (dif2(i)**2.0_dp)
        end do
        delta=(cuad**(0.5_dp))
    end do

    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++    Completa Matriz Real con Ceros    +++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine m_zeros(A,n,m)
    use constants
    implicit none
    integer(dp),intent(in) :: n,m
    real(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j

    do i=1,n,1
        do j=1,m
            A(i,j) = 0.0_dp
        end do
    end do
    
    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++    Completa Vector Real con Ceros    +++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine v_zeros(A,n)
    use constants
    implicit none
    integer(dp),intent(in) :: n
    real(dp),dimension(n),intent(inout) :: A
    integer :: i

    do i=1,n,1
        A(i) = 0.0_dp
    end do
    
    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++    Completa Matriz Compleja con Ceros    +++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine m_zeros_complex(A,n,m)
    use constants
    implicit none
    integer(dp),intent(in) :: n,m
    complex(dp),dimension(n,m),intent(inout) :: A
    integer :: i,j

    do i=1,n,1
        do j=1,m
            A(i,j) = (0.0_dp,0.0_dp)
        end do
    end do
    
    return
end subroutine
!
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++    Completa Vector Complejo con Ceros    +++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine v_zeros_complex(A,n)
    use constants
    implicit none
    integer(dp),intent(in) :: n
    complex(dp),dimension(n),intent(inout) :: A
    integer :: i

    do i=1,n,1
        A(i) = (0.0_dp,0.0_dp)
    end do
    
    return
end subroutine
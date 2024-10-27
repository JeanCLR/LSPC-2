!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++      M A I N   P R O G R A M       ++++++++++++++++
!++++++++++++++++                                    ++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  Solucion Numerica para 'n_r', 'm_max' y 'vector k'
!
!
program main
    use constants
    implicit none
    integer(dp) :: i,j,l
    integer(dp) :: n,m,mp,pa
    real(dp) :: radio, temperatura, funcion_contorno, func_k
    real(dp),dimension(:),allocatable :: k, r, th, sol_aux
    complex(dp),dimension(:,:),allocatable :: coef
    

    write(*,*) '---------------------------------------------'
    write(*,*) '-------------- PROGRAMA LSPC-2 --------------'
    write(*,*) '---------------------------------------------'
    write(*,*) ' '
    write(*,*) 'Ingrese la cantidad de materiales: '
    read(*,*) n
    write(*,*) 'Ingrese el valor de m_max: '
    read(*,*) m

    write(*,*) '¿Ingresar vector de cond. termica? (1=si/0=no)'
    read(*,*) l

    if (l/=0 .AND. l/=1) then
        return
    end if
    
    mp = m+1
    radio = 1.0_dp
    pa = 101_dp
    
    allocate(k(n))
    allocate(r(0:2*n))
    allocate(th(pa))
    allocate(coef(mp,2*n))
    allocate(sol_aux(pa))

    
    call linspace(0.0_dp, radio, r, (2*n)+1)
    call linspace(0.0_dp, 2*pi, th, pa)

    open(20,file='file1.dat',status='replace')
    open(30,file='file2.dat',status='replace')
    open(40,file='file3.dat',status='replace')
    open(50,file='file4.dat',status='replace')
    
    do i=0,2*n,1
        write(20,'(F11.6)') r(i)
        if ((mod(i,2)==0) .AND. (i/=(2*n))) then
            select case (l)
                case (0)
                    k((i/2)+1) = func_k(r(i))
                case (1)
                    write(*,*) 'k(',(i/2)+1,')'
                    read(*,*) k((i/2)+1)
            end select
        end if
    end do

    do i=1,pa,1
        write(30,'(F11.6)') th(i)
    end do

    do i=1,n,1
        write(40,'(F11.6)') k(i)
    end do
    
    call get_coef(n,mp,k,coef)
    
    write(50,*) '#  Matriz de temperaturas T[i,j] ---> T(r_i,th_j)'
    do i=1,2*n,1
        do j=1,pa,1
            sol_aux(j) = temperatura(coef,mp,2*n,i,th(j))
        end do
        !write(50,*) sol_aux
        write(50,'(999(F11.6,X))') (sol_aux(l), l=1,pa)
    end do

    do j=1,pa,1
        sol_aux(j) = funcion_contorno(th(j))
    end do
    write(50,'(999(F11.6,X))') (sol_aux(l), l=1,pa)
    
    deallocate(k)
    deallocate(r)
    deallocate(th)
    deallocate(coef)
    deallocate(sol_aux)

    close(20)
    close(30)
    close(40)
    close(50)

end program main
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
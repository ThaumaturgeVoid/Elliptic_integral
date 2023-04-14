module moduleMELPsubroutine
    implicit none
    Real(8),parameter :: PI=3.141592653589793238D0
    contains
    subroutine melp(k, f, fout,eout)
        real(8),intent(in)  :: k, f
        real(8),intent(out) :: fout,eout
        real(8) :: y, e, ff1, fff1, ff2, fff2
        integer :: n, m, i, j
        If ((k<0.d0) .Or. (k>1.d0)) Then
            stop "k must in 0 <= k <= 1"
            Return
        End If
        ! pi = 3.1415926
        y = dabs(f)
        n = 0
        do while (y>=pi) 
            n = n + 1
            y = y - pi
        End do
        e = 1.d0
        If (y>=pi/2.d0) Then
            n = n + 1
            e = -1.d0
            y = pi - y
        End If
        If (n==0) Then
            call fk(k, y, ff1, ff2)
            ! ff = fk1(k, y)
            fout = ff1
            eout = ff2
            Return
        Else
            call fk(k, pi/2.d0, ff1, ff2)
            ! ff = fk1(k, pi/2.d0)
            call fk(k, y, fff1, fff2)
            ! fff = fk1(k, y)
            fout = 2.d0*n*ff1 + e*fff1
            eout = 2.d0*n*ff2 + e*fff2
            Return
        End If
    End subroutine melp

    subroutine fk(k, y,fk1, fk2)
        real(8) :: t(5), c(5)
        real(8),intent(in) :: k, y
        real(8),intent(out) :: fk1, fk2
        real(8) ::  a, b,  h, s, aa, bb,  x 
        real(8) :: f1, g1, p1, q1, w1
        real(8) :: f2, g2, p2, q2, w2
        integer :: m, i, j
        logical :: fk1continue=.false.,fk2continue=.false.
        t=[-0.9061798459d0, -0.5384693101d0, 0.d0, 0.5384693101d0, 0.9061798459d0]
        c=[0.2369268851d0, 0.4786286705d0, 0.5688888889d0, 0.4786286705d0, 0.2369268851d0]
        a = 0.d0
        b = y
        m = 1
        s = (b-a)*0.02d0
        h = (b-a)/m

        p1 = 0.d0
        p2 = 0.d0

        g1 = 0.d0
        g2 = 0.d0
        Do i = 1, m
            aa = a + (i-1)*h
            bb = a + i*h

            w1 = 0.d0
            w2 = 0.d0

            Do j = 1, 5
                x = ((bb-aa)*t(j)+(bb+aa))/2.d0

                f1 = 1.d0/dsqrt(1.0-k*k*dsin(x)*dsin(x))
                w1 = w1 + f1*c(j)
                f2 = dsqrt(1.0-k*k*dsin(x)*dsin(x))
                w2 = w2 + f2*c(j)

            End Do

            g1 = g1 + w1
            g2 = g2 + w2

        End Do

        g1 = g1*h/2.d0
        q1 = dabs(g1-p1)/(1.0+dabs(g1))
        if((q1>=1.0D-10) .And. (dabs(h)>dabs(s))) fk1continue=.True.

        g2 = g2*h/2.d0
        q2 = dabs(g2-p2)/(1.0+dabs(g2))
        if((q2>=1.0D-10) .And. (dabs(h)>dabs(s))) fk2continue=.True.

        do while ((fk1continue .eqv. .True.) .Or. (fk2continue .eqv. .True.))
            h = (b-a)/m

            g1 = 0.d0
            g2 = 0.d0

            Do i = 1, m
                aa = a + (i-1)*h
                bb = a + i*h

                w1 = 0.d0
                w2 = 0.d0

                Do j = 1, 5
                    x = ((bb-aa)*t(j)+(bb+aa))/2.d0

                    f1 = 1.d0/dsqrt(1.0-k*k*dsin(x)*dsin(x))
                    w1 = w1 + f1*c(j)
                    f2 = dsqrt(1.0-k*k*dsin(x)*dsin(x))
                    w2 = w2 + f2*c(j)

                End Do

                g1 = g1 + w1
                g2 = g2 + w2

            End Do

            g1 = g1*h/2.d0
            q1 = dabs(g1-p1)/(1.d0+dabs(g1))
            if((q1>=1.0D-10) .And. (dabs(h)>dabs(s))) Then
                p1 = g1
            else
                fk1continue=.false.
                fk1 = g1
            end if

            g2 = g2*h/2.d0
            q2 = dabs(g2-p2)/(1.d0+dabs(g2))
            if((q2>=1.0D-10) .And. (dabs(h)>dabs(s))) Then
                p2 = g2
            else
                fk2continue=.false.
                fk2 = g2
            end if
            m = m + 1
        End do
        Return
    End subroutine fk
end module moduleMELPsubroutine
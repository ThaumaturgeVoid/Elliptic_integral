module moduleMELPfunction
    implicit none

    contains
    Function melp1(k, f)
        real(8),intent(in) :: k, f
        real(8) :: melp1, y, ff, fff, e
        integer :: n, m, i, j
        Real(8),parameter :: pi=3.141592653589793238D0
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
            ff = fk1(k, y)
            melp1 = ff
            Return
        Else
            ff = fk1(k, pi/2.d0)
            fff = fk1(k, y)
            melp1 = 2.d0*n*ff + e*fff
            Return
        End If
    End Function melp1

    Function melp2(k, f)
        real(8),intent(in) :: k, f
        real(8) :: melp2, y, ff, fff, e
        integer :: n, m, i, j
        Real(8),parameter :: pi=3.141592653589793238D0
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
            ff = ek1(k, y)
            melp2 = ff
            Return
        Else
            ff = ek1(k, pi/2.d0)
            fff = ek1(k, y)
            melp2 = 2.d0*n*ff + e*fff
            Return
        End If
    End Function melp2

    Function fk1(k, y)
        real(8) :: t(5), c(5)
        real(8),intent(in) :: k, y
        real(8) :: fk1, a, b, f, g, s, p, h, aa, bb, w, x, q
        integer :: n, m, i, j
        t=[-0.9061798459d0, -0.5384693101d0, 0.d0, 0.5384693101d0, 0.9061798459d0]
        c=[0.2369268851d0, 0.4786286705d0, 0.5688888889d0, 0.4786286705d0, 0.2369268851d0]
        a = 0.d0
        b = y
        m = 1
        s = (b-a)*0.02d0
        p = 0.d0
        h = (b-a)/m
        g = 0.d0
        Do i = 1, m
            aa = a + (i-1)*h
            bb = a + i*h
            w = 0.d0
            Do j = 1, 5
                x = ((bb-aa)*t(j)+(bb+aa))/2.d0
                fk1 = 1.d0/dsqrt(1.0-k*k*dsin(x)*dsin(x))
                w = w + fk1*c(j)
            End Do
            g = g + w
        End Do
        g = g*h/2.d0
        q = dabs(g-p)/(1.0+dabs(g))
        do while ((q>=1.0D-10) .And. (dabs(h)>dabs(s)))
            h = (b-a)/m
            g = 0.d0
            Do i = 1, m
                aa = a + (i-1)*h
                bb = a + i*h
                w = 0.d0
                Do j = 1, 5
                    x = ((bb-aa)*t(j)+(bb+aa))/2.d0
                    fk1 = 1.d0/dsqrt(1.d0-k*k*dsin(x)*dsin(x))
                    w = w + fk1*c(j)
                End Do
                g = g + w
            End Do
            g = g*h/2.d0
            q = dabs(g-p)/(1.d0+dabs(g))
            p = g
            m = m + 1
        End do
        fk1 = g
        Return
    End Function fk1

    Function ek1(k, y)
        real(8) :: t(5), c(5)
        real(8),intent(in) :: k, y
        real(8) :: ek1, a, b, f, g, s, p, h, aa, bb, w, x, q
        integer :: n, m, i, j
        t=[-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459]
        c=[0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851]
        a = 0.0
        b = y
        m = 1
        s = (b-a)*0.02
        p = 0.0
        h = (b-a)/m
        g = 0.0
        Do i = 1, m
            aa = a + (i-1)*h
            bb = a + i*h
            w = 0.0
            Do j = 1, 5
                x = ((bb-aa)*t(j)+(bb+aa))/2.0
                ek1 = sqrt(1.0-k*k*sin(x)*sin(x))
                w = w + ek1*c(j)
            End Do
            g = g + w
        End Do
        g = g*h/2.0
        q = abs(g-p)/(1.0+abs(g))
        do while ((q>=1.0D-10) .And. (abs(h)>abs(s)))
            h = (b-a)/m
            g = 0.0
            Do i = 1, m
                aa = a + (i-1)*h
                bb = a + i*h
                w = 0.0
                Do j = 1, 5
                    x = ((bb-aa)*t(j)+(bb+aa))/2.0
                    ek1 = sqrt(1.0-k*k*sin(x)*sin(x))
                    w = w + ek1*c(j)
                End Do
                g = g + w
            End Do
            g = g*h/2.0
            q = abs(g-p)/(1.0+abs(g))
            p = g
            m = m + 1
        End do
        ek1 = g
        Return
    End Function ek1
end module moduleMELPfunction
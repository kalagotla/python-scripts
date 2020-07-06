def pt(n):
    t = [[1]]
    for i in range(1,n):
        r = [1]
        for j in range(1,i):         
            pr = t[i-1]
            r.append(pr[j-1] + pr[j])
        r.append(1)
        t.append(r)    

    me = max(t[j])

    w = len(str(me))

    for i in range(n):
        nf = 2*n-1
        bf = nf//2 - i
        l = " "*bf*w
        for j in range(i+1):
            f = "{{:{}}}".format(w)
            l += f.format(t[i][j]) + " "*w
        print(l)

pt(5)
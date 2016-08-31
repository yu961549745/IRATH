IRATH:=module()
    option  package;
    export  findTanhSolutions;
    local   nomalizeEqn,
            getODE,
            simplifyEqn,
            findDiffOrder,
            findElemOrder,
            findItemOrder,
            findOrders,
            findInflexion,
            solveAllOrder,
            solveForOrder,
            selectSolutions,
            printSolutions,
            solveForPS,
            transSolution,
            subsIeq,
            subsSol,
            checkSolution,
            isEquationOf;

    (*
     * 求解函数
     *)
    findTanhSolutions:=proc(eqi)
        local eq,ms,params,oeq,st;
        st:=time();
        
        printf("输入的方程为");
        eq:=nomalizeEqn(eqi);
        oeq:=eq;
        print(eq);

        # printf("方程中的参数为");
        params:=indets(eq,name) minus {t,x};
        # print(params);
        
        # printf("行波变换后的方程为");
        eq:=getODE(eq);
        # print(eq);

        # printf("方程各项阶数为");
        ms:=findOrders(eq);
        # print(ms);

        # printf("拐点为");
        ms:=findInflexion(ms);
        # print(ms);

        if evalb(ms={0}) then
            printf("无解\n");
            printf("时间已过 %f 秒\n",time()-st);
            return;
        end if;
        
        # printf("开始求解\n");
        solveAllOrder(eq,ms,params,oeq);

        printf("时间已过 %f 秒\n",time()-st);
    end proc:

    (*
     * 标准化方程
     * 方程变量是u(x,t)
    *)
    nomalizeEqn:=proc(ee)
        local e;
        e:=expand(numer(ee));
        return e;
    end proc:

    (*
     * 进行行波变换，将偏微分方程转化为常微分方程
     * 要求函数为u(x,t),替换为u(xi),
     * xi=k*(x-ct)
    *)
    getODE:=proc(eqi)
        local eq;
        eq:=subs(u(x,t)=u(x,t,xi),eqi);
        eq:=PDETools[dsubs](diff(u(x,t,xi),t)=-c*k*diff(u(x,t,xi),xi),
        diff(u(x,t,xi),x)=k*diff(u(x,t,xi),xi),
        eq);
        eq:=eval(subs(u(x,t,xi)=u(xi),eq));
        return simplifyEqn(eq);
    end proc:

    (*
     * 化简方程
    *)
    simplifyEqn:=proc(ee)
        local e,v,p;
        e:=numer(ee);
        e:=factor(e);
        if type(e,`*`) then
            p:=1;
            for v in e do
                if type(v,`+`) then
                    p:=p*v;
                 elif type(v,`^`) and type(expand(op(1,v)),`+`) then
                    p:=p*(simplifyEqn(op(1,v))^op(2,v));
                end if;
            end do;
            p:=expand(p);
        else
            p:=expand(e);
        end if;
        return p;
    end proc:

    (*
     * 计算导数项阶数
    *)
    findDiffOrder:=proc(ee)
        local e;
        e:=convert(ee,'D');
        e:=op(0,op(0,e));
        if evalb(op(0,e)=`@@`) then
            op(2,e);
        else
            1;
        end if;
    end proc:

    (*
     * 计算元素的阶数
    *)
    findElemOrder:=proc(ee)
        local e;
        if evalb(op(0,ee)='u') then
            e:=[1,0];
        elif evalb(op(0,ee)='diff') then
            e:=[1,findDiffOrder(ee)];
        elif type(ee,`^`) then
            e:=findElemOrder(op(1,ee))*op(2,ee);
        else
            e:=[0,0];
        end if;
        return e;
    end proc:

    (*
     * 计算方程某一项的阶数
    *)
    findItemOrder:=proc(ee)
        local e,p;
        if type(ee,`*`) then
            p:=[0,0];
            for e in ee do
                p:=p+findElemOrder(e);
            end do;
        else
            p:=findElemOrder(ee);
        end if;
        return p;
    end proc:

    (*
     * 计算方程每一项的阶数
     * 假设至少有两项
     * 假设指数上不含参数
    *)
    findOrders:=proc(ee)
        findItemOrder~([op(ee)]);
    end proc:

    (*
     * 计算拐点
    *)
    findInflexion:=proc(mss::list(list(integer)))
        local ms,i,j,n,mm,m,fun;
        fun:=(x,m)->x[1]*m+x[2];
        mm:={};
        ms:={op(mss)};
        n:=numelems(ms);
        for i from 1 to n-1 do
            for j from i+1 to n do
                if evalb(ms[i][1]<>ms[j][1]) then
                    m:=-(ms[i][2]-ms[j][2])/(ms[i][1]-ms[j][1]);
                    if evalb(fun(ms[i],m)=max(map(fun,ms,m))) then
                        mm:=mm union {m};
                    end if;
                end if;
            end do;
        end do;
        return mm;
    end proc:

    (*
     * 对所有阶数进行求解
    *)
    solveAllOrder:=proc(eqi,mms::set(rational),params::set(name),oeq)
        local trans,tr,eq,ms,m;
        # 计算全部变换
        trans:=map(x->sign(x)/denom(x),mms);
        for tr in trans do
            eq:=eval(subs(u(xi)=u(xi)^tr,eqi));
            eq:=simplifyEqn(eq);
            ms:=findOrders(eq);
            m:=select(type,findInflexion(ms),posint);
            m:=max(m);
            printf("m=%a,方程为:",tr*m);
            print(eq);
            if not type(ms,list(list(integer))) then
                printf("方程不是多项式，无解\n");
                next;
            end if;
            solveForOrder(eq,m,params,tr,oeq);
        end do;
    end proc:

    (*
     * 对特定阶数进行求解
    *)
    solveForOrder:=proc(eqi,m::posint,params::set(name),tr::rational,oeq)
        local eq,f,PS,vars,sols;
        f:=add(seq(a[i]*tanh(xi)^i,i=0..m));
        vars:=[seq(a[i],i=0..m)];
        eq:=eval(subs(u(xi)=f,eqi));
        eq:=subs(tanh(xi)=T,eq);
        eq:=collect(eq,T);
        PS:=PolynomialTools[CoefficientList](eq,T,termorder=reverse);
        sols:=solveForPS(PS,m,params,vars);
        printSolutions(sols,m,tr);
    end proc:

    (*
     * 筛选解
     * 删去k,c=0的解
     * 删除参数为0的解
     * 删除a[1]..a[m]都为0的解
    *)
    selectSolutions:=proc(solss::{set,list},
                         params::set(name),m::posint)
        local zs,ps,sols,as;
        as:={seq(a[i],i=1..m)};
        sols:=convert(solss,set);
        sols:=select(type,sols,equation);
        ps:=params union {k,c};
        zs:=lhs~(select(_x->type(_x,`=`(name,0)),sols));
        if evalb( (zs intersect ps) <> {} ) then
            return false;
        else
            return not evalb(as subset zs);
        end if; 
    end proc:

    (*
     * 输出解
    *)
    printSolutions:=proc(sols::set({list,set}),
                        m::posint,tr::rational)
        local sol;
        if evalb(sols={}) then
            printf("无解\n");
            return;
        end if;
        printf("共有%d个解",numelems(sols));
        print(u(xi)=add(seq(a[i]*tanh(xi)^i,i=0..m))^tr,xi=k*(x-ct)+xi[0]);
        printf("其中");
        for sol in sols do
            print(expand(simplify(sol)));
        end do;
    end proc:

    (*
     * 求解PS
     * 基于csolve的版本
     * 只有剩余方程不能用csolve求解时，才使用solve求解
    *)
    solveForPS:=proc(PS,m::posint,params::set(name),vars::list(name))
        local pa,pp,sola,sol,vs,res,np,nsol,ssol,nres,sssol;
        res:={};# 所有解
        pa:=PS[1..(m+1)];
        pp:=PS[(m+1)..-1];
        sola:=[csolve(pa,{vars[],k,c})];# 求解前m+1个方程
        sola:=select(selectSolutions,sola,params,m);# 除去平凡的解
        for sol in sola do
            # 代入化简
            np:=simplify(subs(op(sol),pp));
            np:=simplifyEqn~(np);# 化简方程，去除非零项
            np:=remove(type,np,0);
            
            # 如果剩余的方程全0，则不需要继续求解
            if evalb(np=[]) then
                res:=res union {sol};
                next;
            end if;
            
            # 求解剩余方程
            try
                nsol:=[csolve(np,indets(np,name))];
            catch:
                nsol:=[RealDomain[solve](np,indets(np,name))];
            end try;
            
            nsol:=select(selectSolutions,nsol,params,m);
            
            for ssol in nsol do
                # 合并解，因为前面的解，在后面已经被带入
                # 所以，后面的解不会包含前面已有的符号
                # 因此，直接合并即可
                try # 有可能出现分母为0的情况
                    sssol:=simplify(subs(op(ssol),sol));
                    res:=res union {sssol union ssol};
                catch:
                    next;
                end try;
            end do;
        end do;
        
        res:=select(selectSolutions,res,params,m);
        
        # printf("未转换的解");
        # print~(res);
        
        nres:={};
        map(transSolution,res,'nres',params,vars);
        res:=select(selectSolutions,nres,params,m);
        
        return res;
    end proc:



    (*
     * 转化csolve的结果
    *)
    transSolution:=proc(res::set(equation),nres::evaln(set),
                       params::set(name),vars::list(name))
        local r,rv,rp,np,rrv,rrp,sp,newSol;
        # 先求解a[0]..a[m],k,c
        sp:=select(isEquationOf,res,{vars[],k,c});
        sp:=map(_x->simplifyEqn(lhs(_x)-rhs(_x)),sp);# 化简方程，去除非零项
        rv:=RealDomain[solve](sp,[vars[],k,c]);

        # 如果求解失败则直接原样返回
        if evalb(rv=[]) then
            nres:=eval(nres) union {[res[]]};
            WARNING("solve求解失败");
            return;
        end if;

        # 再求解参数
        for rrv in rv do
            try
                np:=simplify(subs(op(rrv),res));
            catch:
                next;
            end try;
            np:=map(_x->simplifyEqn(lhs(_x)-rhs(_x)),np);# 化简方程，去除非零项
            np:=select(isEquationOf,np,params);
            rp:=RealDomain[solve](np,[params[]]);
            for rrp in rp do
                rrp:=subsIeq~(rrp);
                try
                    newSol:=[op(subsSol(rrv,rrp)),rrp[]];
                    newSol:=simplify(rationalize(newSol));
                    nres:=eval(nres) union {newSol};
                catch:
                    next;
                end try;
            end do;
        end do;
        return;
    end proc:

    (*
     * 消去不等式范围约束
    *)
    subsIeq:=proc(ee)
        local v;
        if type(ee,equation) then
            return ee;
        else
            v:=indets(ee);
            if evalb(nops(v)=1) then
                v:=v[1];
                return (v=v);
            else
                return ee;
            end if;
        end if;
    end proc:

    (*
     * 解集带入
     * sol2代入sol1
    *)
    subsSol:=proc(sol1,sol2)
        local sol;
        sol:=select(type,sol2,equation);
        expand(subs(op(sol),sol1));
    end proc:

    (*
     * 代回原方程进行验证
    *)
    checkSolution:=proc(ssol,oeq,m,tr)
        local uxi,sol,r;
        sol:=select(type,ssol,equation);
        r:=subs(u(x,t)=add(seq(a[i]*tanh(xi)^i,i=0..m))^tr,oeq);
        r:=subs(xi=k*(x-c*t),r);
        r:=subs(op(sol),r);
        r:=simplify(r);
        evalb(r=0);
    end proc:

    (*
     * 选择关于给定参数的方程
    *)
    isEquationOf:=proc(eq,vars::set)
        local vs;
        vs:=indets(eq,name);
        evalb( (vs intersect vars) <> {} );
    end proc:

end module:
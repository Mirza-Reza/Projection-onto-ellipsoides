using LinearAlgebra
function Quadratic(x,A,b)
    return dot(x,A*x)+2*dot(b,x)
end


 function GenAlg(x₀,A,b,α,a ;ϵ=1e-6, max_itr=999 )
    xₖ=x₀
    gₐ=A*a
    qa_α = Quadratic(a,A,b)-α
    uₖ =A*xₖ+b
    vₖ =a-xₖ 
    uₖdotvₖ=dot(uₖ,vₖ)
    vₖdotvₖ=dot(vₖ,vₖ)
    normuₖ=norm(uₖ)
    normvₖ=sqrt(vₖdotvₖ)
    Mₖ =uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
    E=1-Mₖ 
    k=0
    while k<= max_itr && E>=ϵ
        zₖ =uₖ-((uₖdotvₖ)/(vₖdotvₖ))*vₖ
        zₖdotzₖ=dot(zₖ,zₖ)
        Avₖ=A*vₖ 
        Azₖ=A*zₖ 
        vₖdotAvₖ=dot(vₖ,Avₖ)
        zₖdotAzₖ=dot(zₖ,Azₖ)
        vₖdotAzₖ=dot(vₖ,Azₖ)
        M1=(vₖdotAvₖ)/(vₖdotvₖ)
        M2=(zₖdotAzₖ)/(zₖdotzₖ)
        M3=((M1-M2)^2 )
        M4=4*(vₖdotAzₖ)^2/(vₖdotvₖ*zₖdotzₖ)
        M5=sqrt(M3+M4)
        ρₖ=[M1+M2+M5]/2 
        γₖ=1/ρₖ
        cₖ=xₖ-( γₖ.*uₖ)  #the center of maximal 2-d inside ball
        #  Step3  ==== Calculation of  xk+1  ==============================
        wₖ=cₖ-a
        # ηk =============================
        Awₖ=A*wₖ
        wₖdotAwₖ= dot(wₖ,Awₖ)

         
        N1=dot(gₐ,wₖ)/(wₖdotAwₖ)
        N2=(N1).^2
        N3=(qa_α)./(wₖdotAwₖ)
        ηₖ=-N1-sqrt(N2-N3)  #η is a number too close to 1.
        xₖ=a+(ηₖ.*wₖ)  
        uₖ =A*xₖ+b
        vₖ =a-xₖ 
        uₖdotvₖ=dot(uₖ,vₖ)
        vₖdotvₖ=dot(vₖ,vₖ)
        normuₖ=norm(uₖ)
        normvₖ=sqrt(vₖdotvₖ)
        Mₖ =uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
        E=1-Mₖ 
        k+=1
        #println("It :$k,  x=$xk" )



    end
    return xₖ,k,E
end
     
     N=10
     i = 1:1:N
     i² =i.^2
     a =(i².*10).+1
     A=Diagonal(i²)
     b=zeros(N)
     α=385
     x₀=sqrt(α/Quadratic(a,A,b)).*a 
    
     xₖ,k,E=  GenAlg(x₀,A,b,α,a,ϵ=1e-6)
# Tested on example 1 i section3 and with 111 iteration gives exactlly same results as Dia's paper.
using LinearAlgebra
using PolynomialRoots
function Quadratic(x,A,b)
    return dot(x,A*x)+2*dot(b,x)
end

function  pickupRoot(β,h₁,h₂,λ₁,λ₂)   
    p = [0.,0]
    #solve Quadratic in lemma 1 to find p₁ then p₂ with following formula
    #P=((λ₁*(λ₁-λ₂)^2))*(p₁)^4+(2*λ₁*λ₂*h₁*(λ₁-λ₂))*(p₁)^3+((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))*(p₁^2)-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)*(p₁)-(β*λ₂^2*h₁^2)
    a4=((λ₁*(λ₁-λ₂)^2))
    a3=(2*λ₁*λ₂*h₁*(λ₁-λ₂))
    a2=((λ₁*λ₂^2*h₁^2)+(λ₁^2*λ₂*h₂^2)-(β*(λ₁-λ₂)^2))
    a1=-2*β*(λ₁*λ₂*h₁-λ₂^2*h₁)
    a0=-(β*λ₂^2*h₁^2)
    S=roots([a0, a1, a2, a3, a4])
    #choose just Real roots =============================
    indexreal=findall(x->abs.(x)<1e-12,imag.(S))
    Sol=real.(S[indexreal])
 
    #choose the right root ============================================
    !isempty(Sol) ? nothing :  @error "No solution for the quartic equation"
    for  p₁ in Sol
        sol1 = p₁
        μ₁=(h₁-p₁)/(λ₁*p₁)
        if μ₁<=0.
            continue
        else
            p₂=(λ₁*h₂)/((λ₁-λ₂)*p₁+λ₂*h₁)*p₁
            μ₂=(h₂-p₂)/(λ₂*p₂)
            if  μ₂<=0. && !(μ₂ ≈ μ₁) 
                continue
            end
        end
        sol2 = p₂
        p = [sol1,sol2]
    end
    return   p
end

function compute_Alg1_γk(vₖdotAvₖ,vₖdotvₖ,vₖdotAzₖ,zₖdotAzₖ,zₖdotzₖ)
    M1=(vₖdotAvₖ)/(vₖdotvₖ)
    M2=(zₖdotAzₖ)/(zₖdotzₖ)
    M3=((M1-M2)^2 )
    M4=4*(vₖdotAzₖ)^2/(vₖdotvₖ*zₖdotzₖ)
    M5=sqrt(M3+M4)
    ρₖ=0.5*(M1+M2+M5)
    γₖ=1/ρₖ
    return γₖ
end

function compute_Alg1_ηk(gₐ,qa_α,wₖ,wₖdotAwₖ)
    N1=dot(gₐ,wₖ)/(wₖdotAwₖ)
    N2=(N1).^2
    N3=(qa_α)./(wₖdotAwₖ)
    ηₖ=-N1-sqrt(N2-N3)
    return ηₖ
end





 function GenAlg(x₀,A,b,α,a ;ϵ=1e-6, max_itr=200)
    xₖ=x₀
    gₐ=A*a + b
    qa_α = Quadratic(a,A,b)-α
    uₖ =A*xₖ+b
    vₖ =a-xₖ 
    uₖdotvₖ=dot(uₖ,vₖ)
    vₖdotvₖ=dot(vₖ,vₖ)
    normuₖ=norm(uₖ)
    normvₖ=sqrt(vₖdotvₖ)
    Mₖ=uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
    E=1-Mₖ 
    k=0
    called_alg = :Alg2
    while k<= max_itr && E>=ϵ

        zₖ =uₖ-((uₖdotvₖ)/(vₖdotvₖ))*vₖ
        normzₖ=norm(zₖ)
        zₖdotzₖ=dot(zₖ,zₖ)
        uₖdotzₖ=dot(uₖ,zₖ)
        Avₖ=A*vₖ 
        Azₖ=A*zₖ 
        vₖdotAvₖ=dot(vₖ,Avₖ)
        zₖdotAzₖ=dot(zₖ,Azₖ)
        vₖdotAzₖ=dot(vₖ,Azₖ)
        zₖdotAvₖ=dot(zₖ,Avₖ)
       #= 
        cₖ=xₖ-( γₖ⁽¹⁾.*uₖ)  #the center of maximal 2-d inside ball
        #  Step3  ==== Calculation of  xk+1  ==============================
        wₖ=cₖ-a
        wₖ=xₖ-( γₖ⁽¹⁾.*uₖ)-a
        # ηk =============================
        Awₖ=A*wₖ
        wₖdotAwₖ= dot(wₖ,Awₖ)
         
        N1=dot(gₐ,wₖ)/(wₖdotAwₖ)
        N2=(N1).^2
        N3=(qa_α)./(wₖdotAwₖ)
        ηₖ=-N1-sqrt(N2-N3)  #η is a number too close to 1.
        xₖ=a+(ηₖ.*wₖ)  =#

       #Alg2================================================
       aₗ=[normvₖ,0]
       Aₖ=[(vₖdotAvₖ)/(normvₖ^2)  (vₖdotAzₖ)/(normvₖ*normzₖ) ; 
             (zₖdotAvₖ)/(normvₖ*normzₖ)  (zₖdotAzₖ)/(normzₖ^2)]
       bₖ=[(uₖdotvₖ)/(normvₖ);(uₖdotzₖ)/(normzₖ)]
       D, Q = eigen(Aₖ)
       Qbₖ=Q*bₖ
       D⁻¹Qbₖ=D.\Qbₖ
       #Aₖ=Q'*D*Q
       β=dot(Qbₖ,D⁻¹Qbₖ)
       aₕ=Q*aₗ+D⁻¹Qbₖ
       #Solve Quadratic=========================================
       h₁=aₕ[1]
       h₂=aₕ[2]
       λ₁=D[1]
       λ₂=D[2]
       p  = pickupRoot(β,h₁,h₂,λ₁,λ₂)
        # println(p)
        #println(" p₁=$p₁,p₂=$p₂")
        #return   p₁,p₂
    aₕ⁽ᵖ⁾= p       
    aₗ⁽ᵖ⁾=Q'*(aₕ⁽ᵖ⁾- D⁻¹Qbₖ)
      # xₖ₊₁ ===========================================================#
    ηₖ⁽²⁾ = 1 - (aₗ⁽ᵖ⁾[1]/normvₖ- (uₖdotvₖ/vₖdotvₖ)*aₗ⁽ᵖ⁾[2]/normzₖ)
    γₖ⁽²⁾=(- aₗ⁽ᵖ⁾[2]/normzₖ)/(ηₖ⁽²⁾)
    #    @show γₖ⁽²⁾
        # γₖ⁽²⁾= -1
       if γₖ⁽²⁾ < 0. || abs(γₖ⁽²⁾) <= ϵ
            γₖ⁽²⁾ = compute_Alg1_γk(vₖdotAvₖ,vₖdotvₖ,vₖdotAzₖ,zₖdotAzₖ,zₖdotzₖ)
            # @show γₖ⁽²⁾
            ωₖ⁽²⁾=-(γₖ⁽²⁾*uₖ)-vₖ
            wₖdotAwₖ = dot(ωₖ⁽²⁾,A*ωₖ⁽²⁾)
            ηₖ⁽²⁾ = compute_Alg1_ηk(gₐ,qa_α,ωₖ⁽²⁾,wₖdotAwₖ)
       else
            ωₖ⁽²⁾=-1*(γₖ⁽²⁾*uₖ)-vₖ
       end
       xₖ=a.+(ηₖ⁽²⁾*ωₖ⁽²⁾)
       uₖ =A*xₖ+b
       vₖ =a-xₖ 
       uₖdotvₖ=dot(uₖ,vₖ)
       vₖdotvₖ=dot(vₖ,vₖ)
       normuₖ=norm(uₖ)
       normvₖ=sqrt(vₖdotvₖ)
       Mₖ=uₖdotvₖ/(normuₖ* normvₖ)  #Mk is qoutient in tolerance condition
       E=1-Mₖ 
       k+=1
    end
    
    return xₖ,k,E
end
    
#=
     N=10
     i = 1:1:N
     i² =i.^2
     a =(i².*10).+1
     A=Diagonal(i²)
     b=zeros(N)
     α=385
     x₀=sqrt(α/Quadratic(a,A,b)).*a 
      xₖ,k,E=  GenAlg(x₀,A,b,α,a,ϵ=1e-6)
=#
     N=3
     i = [1,10,100]
     i² =i.^2
     #a =(i².*1).+1
     a=[2,1.01,1.0001]
     A=Diagonal(i²)
     b=zeros(N)
     α=1.0101
      x₀=[0.01,0.01,0.01]
     #x₀=sqrt(α/Quadratic(a,A,b)).*a 
      xₖ,k,E=  GenAlg(x₀,A,b,α,a,ϵ=1e-6)


     #=Alg3
     γ₃ₖ₊ᵢ=γ₃ₖ₊ᵢ⁽²⁾     for i=0,1
     γ₃ₖ₊₂=γ₃ₖ₊ᵢ⁽¹⁾ 
     =#
     #=================================================
max_itr=5
        k=1
       #i=0 
        while k<=max_itr
       
          for i=0:2
                if i<2,
,
,
,
,
            end
         k+=1
        continue  
        end
      =#

      
      #=Alg4============================================================
     m1=1
m2=1
c1=0.1
c2=0.8
max_itr=20
k=1
#i=0 
while k<=max_itr
     b=mod(k,c1+c2) # Remainder after division (modulo operation)
        if b<c1
            println("use Alg 2 to find γₖ⁽²⁾")
        else 
            println(" γₖ=c1*γₖ⁽¹⁾+ c2*γₖ⁽²⁾")
        end
 k+=1
continue  
end
     =#
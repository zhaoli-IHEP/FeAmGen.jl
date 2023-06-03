using Dates, FeAmGen, SymEngine, FeynUtils, Test

@info "basicTest starts @ $(now())"

@testset "vectorized_tensor_product_String" begin
  target_in = Array{String}[ String["a","b","c"], String["d","e","f"] ]
  target_out = String["a,d","a,e","a,f","b,d","b,e","b,f","c,d","c,e","c,f"]
  @test (isempty∘setdiff)( FeAmGen.vectorized_tensor_product_String( target_in... ), target_out )
  target_in = Array{String}[ String["a","b"], String["c","d"], String["e","f"] ]
  target_out = String["a,c,e","a,c,f","a,d,e","a,d,f","b,c,e","b,c,f","b,d,e","b,d,f"]
  @test (isempty∘setdiff)( FeAmGen.vectorized_tensor_product_String( target_in... ), target_out )
end # @testset


@testset "expand_parton" begin
  @test FeAmGen.expand_parton( String["u","parton"], String["d","g"], String["u","d","g"] ) == String["u,u,d,g", "u,d,d,g", "u,g,d,g"]
end # @testset


@testset "sort_proc_str" begin
  @test FeAmGen.sort_proc_str( "u,d,g,d,u", 2 ) == "d,u,d,g,u"
end # @testset


@testset "get_list_quoted_str" begin
  @test FeAmGen.get_list_quoted_str( String["a","b","c"] ) == "[ \"a\",\"b\",\"c\" ]"
end # @testset





@testset "Test generate_kin_relation_v2" begin

  @vars K1, K2, K3, K4, m1, m2, m3, m4, shat

  nn              =   50
  mom_list_str    =   join(["K$ii" for ii ∈ 1:nn], ", ")
  mass2_list_str  =   join(["m$ii" for ii ∈ 1:nn], ", ")
  ver_list_str    =   join(["ver$ii" for ii ∈ 1:(3 * nn - 10)], ", ")
  mom_list        =   [Basic("K$ii") for ii ∈ 1:nn]
  mass2_list      =   [Basic("m$ii") for ii ∈ 1:nn]
  # alias
  Ki = mom_list
  mi = mass2_list

  (eval ∘ Meta.parse)( "@vars $mom_list_str, $mass2_list_str, $ver_list_str")
  half = Basic(1//2)


  #---------------------
  # 1 -> 1
  n_inc = 1
  n_out = 1
  n_tot = n_inc + n_out

  kin_relation = FeAmGen.generate_kin_relation_v2(
    n_inc, n_out,
    mom_list[1:n_tot],
    mass2_list[1:n_tot]
  )
  kin_relation_bench = Dict{Basic, Basic}(
    make_SP(K1,K2) => m1
  )
  for ii ∈ 1:n_tot
    Kii = Ki[ii]
    mii = mi[ii]
    key = make_SP(Kii,Kii)
    kin_relation_bench[key] = mii
  end # for ii
  @test kin_relation == kin_relation_bench
  # 1 -> 1
  #---------------------


  #---------------------
  # 1 -> 2
  n_inc = 1
  n_out = 2
  n_tot = n_inc + n_out

  kin_relation = FeAmGen.generate_kin_relation_v2(
    n_inc, n_out,
    mom_list[1:n_tot],
    mass2_list[1:n_tot]
  )
  kin_relation_bench = Dict{Basic, Basic}(
    # Basic("SP(K1, K2)") => Basic("(1/2)*m1 + (1/2)*m2 + (-1/2)*m3"),
    # Basic("SP(K1, K3)") => Basic("(1/2)*m1 + (-1/2)*m2 + (1/2)*m3"),
    # Basic("SP(K2, K3)") => Basic("(1/2)*m1 + (-1/2)*m2 + (-1/2)*m3"),
    make_SP(K1, K2) => half * m1 + half * m2 - half * m3,
    make_SP(K1, K3) => half * m1 - half * m2 + half * m3,
    make_SP(K2, K3) => half * m1 - half * m2 - half * m3
  )
  for ii ∈ 1:n_tot
    Kii = Ki[ii]
    mii = mi[ii]
    key = make_SP(Kii,Kii)
    kin_relation_bench[key] = mii
  end # for ii
  @test kin_relation == kin_relation_bench
  # 1 -> 2
  #---------------------


  #---------------------
  # 2 -> 1
  n_inc = 2
  n_out = 1
  n_tot = n_inc + n_out

  kin_relation = FeAmGen.generate_kin_relation_v2(
    n_inc, n_out,
    mom_list[1:n_tot],
    mass2_list[1:n_tot]
  )
  kin_relation_bench = Dict{Basic, Basic}(
    # Basic("SP(K1, K2)") => Basic("(-1/2)*m1 + (-1/2)*m2 + (1/2)*m3"),
    # Basic("SP(K1, K3)") => Basic("(1/2)*m1 + (-1/2)*m2 + (1/2)*m3"),
    # Basic("SP(K2, K3)") => Basic("(-1/2)*m1 + (1/2)*m2 + (1/2)*m3")
    make_SP(K1, K2) => - half * m1 - half * m2 + half * m3,
    make_SP(K1, K3) => half * m1 - half * m2 + half * m3,
    make_SP(K2, K3) => - half * m1 + half * m2 + half * m3
  )
  for ii ∈ 1:n_tot
    Kii = Ki[ii]
    mii = mi[ii]
    key = make_SP(Kii,Kii)
    kin_relation_bench[key] = mii
  end # for ii
  @test kin_relation == kin_relation_bench
  # 2 -> 1
  #---------------------


  #---------------------
  # 2 -> 2
  n_inc = 2
  n_out = 2
  n_tot = n_inc + n_out

  kin_relation = FeAmGen.generate_kin_relation_v2(
    n_inc, n_out,
    mom_list[1:n_tot],
    mass2_list[1:n_tot]
  )
  kin_relation_bench = Dict{Basic, Basic}(
    # Basic("SP(K1, K2)") => Basic("(1/2)*(-m1 - m2 + shat)"),
    # Basic("SP(K1, K3)") => Basic("(-1/2)*(-m1 - m3 + ver1)"),
    # Basic("SP(K1, K4)") => Basic("(-1/2)*m2 + (-1/2)*m3 + (1/2)*shat + (1/2)*ver1"),
    # Basic("SP(K2, K3)") => Basic("(-1/2)*m1 + (-1/2)*m4 + (1/2)*shat + (1/2)*ver1"),
    # Basic("SP(K2, K4)") => Basic("(1/2)*m2 + (1/2)*m4 + (-1/2)*ver1"),
    # Basic("SP(K3, K4)") => Basic("(-1/2)*m3 + (-1/2)*m4 + (1/2)*shat")
    make_SP(K1,K2) => half * (- m1 - m2 + shat),
    make_SP(K1,K3) => - half * (- m1 - m3 + ver1),
    make_SP(K1,K4) => - half * m2 - half * m3 + half * shat + half * ver1,
    make_SP(K2,K3) => - half * m1 - half * m4 + half * shat + half * ver1,
    make_SP(K2,K4) => half * m2 + half * m4 - half * ver1,
    make_SP(K3,K4) => - half * m3 - half * m4 + half * shat
  )
  for ii ∈ 1:n_tot
    Kii = Ki[ii]
    mii = mi[ii]
    key = make_SP(Kii,Kii)
    kin_relation_bench[key] = mii
  end # for ii
  @test kin_relation == kin_relation_bench
  # 2 -> 2
  #---------------------


  #---------------------
  # 2 -> 3
  n_inc = 2
  n_out = 3
  n_tot = n_inc + n_out

  kin_relation = FeAmGen.generate_kin_relation_v2(
      n_inc, n_out,
      mom_list[1:n_tot],
      mass2_list[1:n_tot]
  )
  kin_relation_bench = Dict{Basic, Basic}(
    # Basic("SP(K1, K2)") => Basic("(1/2)*(-m1 - m2 + shat)"),
    # Basic("SP(K1, K3)") => Basic("(-1/2)*(-m1 - m3 + ver1)"),
    # Basic("SP(K1, K4)") => Basic("(-1/2)*(-m1 - m4 + ver2)"),
    # Basic("SP(K1, K5)") => Basic("(-1/2)*m3 + (-1/2)*m4 + (1/2)*ver1 + (1/2)*ver2 + SP(K1, K2)"),
    # Basic("SP(K2, K3)") => Basic("(-1/2)*(-m2 - m3 + ver3)"),
    # Basic("SP(K2, K4)") => Basic("(-1/2)*(-m2 - m4 + ver4)"),
    # Basic("SP(K2, K5)") => Basic("(-1/2)*m3 + (-1/2)*m4 + (1/2)*ver3 + (1/2)*ver4 + SP(K1, K2)"),
    # Basic("SP(K3, K4)") => Basic("(1/2)*m1 + (1/2)*m2 + (1/2)*m3 + (1/2)*m4 + (1/2)*m5 + (-1/2)*ver1 + (-1/2)*ver2 + (-1/2)*ver3 + (-1/2)*ver4 - SP(K1, K2)"),
    # Basic("SP(K3, K5)") => Basic("(-1/2)*m3 + (-1/2)*m4 + (-1/2)*m5 + (1/2)*ver2 + (1/2)*ver4 + SP(K1, K2)"),
    # Basic("SP(K4, K5)") => Basic("(-1/2)*m3 + (-1/2)*m4 + (-1/2)*m5 + (1/2)*ver1 + (1/2)*ver3 + SP(K1, K2)")
    make_SP(K1,K2) => half*(- m1 - m2 + shat),
    make_SP(K1,K3) => - half * (- m1 - m3 + ver1),
    make_SP(K1,K4) => - half * (- m1 - m4 + ver2),
    make_SP(K1,K5) => - half*m1 - half*m2 - half*m3 - half*m4 + half*ver1 + half*ver2 + half*shat,
    make_SP(K2,K3) => - half * (- m2 - m3 + ver3),
    make_SP(K2,K4) => - half * (- m2 - m4 + ver4),
    make_SP(K2,K5) => - half*m1 - half*m2 - half*m3 - half*m4 + half*ver3 + half*ver4 + half*shat,
    make_SP(K3,K4) => m1 + m2 + half*m3 + half*m4 + half*m5 - half*ver1 - half*ver2 - half*ver3 - half*ver4 - half*shat,
    make_SP(K3,K5) => - half*m1 - half*m2 - half*m3 - half*m4 - half*m5 + half*ver2 + half*ver4 + half*shat,
    make_SP(K4,K5) => - half*m1 - half*m2 - half*m3 - half*m4 - half*m5 + half*ver1 + half*ver3 + half*shat 
  )
  for ii ∈ 1:n_tot
    Kii = Ki[ii]
    mii = mi[ii]
    key = make_SP(Kii,Kii)
    kin_relation_bench[key] = mii
  end # for ii

  @assert length(kin_relation) == length(kin_relation_bench)
  for key ∈ keys(kin_relation) 
    @test kin_relation[key] == kin_relation_bench[key]
  end # for key
  # 2 -> 3
  #---------------------


  #---------------------
  # 2 -> (n - 2)
  n_inc = 2
  for n_tot ∈ 6:nn   
    n_out = n_tot - n_inc
    sign_list = [ii ≤ n_inc ? (- 1) : 1 for ii ∈ 1:n_tot]

    kin_relation_bench = Dict{Basic,Basic}([ make_SP(Ki[ii],Ki[ii]) => mi[ii] for ii ∈ 1:n_tot ])
    kin_relation_bench[make_SP(K1,K2)] = half*(-m1 - m2 + shat)
    for ii ∈ 3:(n_tot - 1)
      Kii = Ki[ii]
      signii = sign_list[ii]
      key = make_SP(K1,Kii)
      kin_relation_bench[key] = first(sign_list)*signii*half*Basic("-m1-m$ii+ver$(ii-2)")
    end # for ii
    ver_index = n_tot - 2
    for ii ∈ 2:(n_tot - 3)
      for jj ∈ (ii + 1):(n_tot - 1)
        Kii = Ki[ii]
        Kjj = Ki[jj]
        signii = sign_list[ii]
        signjj = sign_list[jj]
        key = make_SP(Kii,Kjj)
        kin_relation_bench[key] = signii*signjj*half*Basic("-m$ii-m$jj+ver$(ver_index)")
        ver_index += 1
      end # for jj
    end # for ii
    @assert ver_index == div( n_tot * (n_tot - 3), 2 )

    for ii ∈ 1:(n_tot - 3)
      Kii = Ki[ii]
      mii = mi[ii]
      signii = sign_list[ii]

      sp_list = Vector{Basic}()
      for jj ∈ 1:(n_tot-1)
        Kjj = Ki[jj]
        signjj = sign_list[jj]
        push!( sp_list, signii*signjj * ( ii == jj ? 0 : make_SP(Kii,Kjj) ) )
      end # for jj
      tmp = -mii-sum(sp_list)

      Knn = Ki[n_tot]
      signnn = sign_list[n_tot]
      key = make_SP(Kii,Knn)
      kin_relation_bench[key] = signii*signnn*(expand∘subs)(tmp,kin_relation_bench)
    end # for ii

    tmp_list = Vector{Basic}()
    for ii ∈ (n_tot-2):n_tot
      Kii = Ki[ii]
      mii = mi[ii]
      signii = sign_list[ii]

      sp_list = Vector{Basic}()
      for jj ∈ 1:(n_tot-3)
        Kjj = Ki[jj]
        signjj = sign_list[jj]
        push!( sp_list, signjj*make_SP(Kjj,Kii) )
      end # for jj

      push!( tmp_list, -mii-signii*sum(sp_list) )
    end # for ii
    tmp = inv( Basic[ 1 1 0; 1 0 1; 0 1 1 ] ) * tmp_list
    tmp = (expand∘subs).(tmp, Ref(kin_relation_bench))
    push!( kin_relation_bench, make_SP(Ki[n_tot-2],Ki[n_tot-1]) => tmp[1] )
    push!( kin_relation_bench, make_SP(Ki[n_tot-2],Ki[n_tot]) => tmp[2] )
    push!( kin_relation_bench, make_SP(Ki[n_tot-1],Ki[n_tot]) => tmp[3] )
    
    kin_relation = FeAmGen.generate_kin_relation_v2(
      n_inc, n_out,
      mom_list[1:n_tot],
      mass2_list[1:n_tot]
    )
    
    @test length(kin_relation) == length(kin_relation_bench)
    for ii ∈ 1:n_tot
      Kii = Ki[ii]
      for jj ∈ (ii+1):n_tot
        Kjj = Ki[jj]
        key = make_SP(Kii,Kjj)
        @test (iszero∘expand)(kin_relation[key]-kin_relation_bench[key])
      end # for jj
    end # for ii

  end # for n_tot
end # @testset



















@info "basicTest ends @ $(now())"


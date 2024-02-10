


########################################################################################
# \begin{quote}
# \hypertarget{header-n3}{%
# \section{\texorpdfstring{Kinematic variables for (\(n-2\))-body final
# state}{Kinematic variables for (n-2)-body final state}}\label{header-n3}}
# \end{quote}
# 
# Considering \(2\to (n-2)\) process, we set up momenta
# \(k_1,k_2,k_3,\dots,k_n\).
# 
# \textbf{According to RAJENDRA KUMAR, Phys.Rev. D2 (1970) 1902-1914}
# 
# For (\(n-2\))-body final state, there Eq. (2.6) shows the \((3(n-2)-4)\)
# independent variables are
# 
# \[s_{i} = \left( k_1 -\sum_{j=3}^{i+2} k_j \right)^2, \quad \text{for}~~i=1,\dots,n-3\label{eq1}\]
# 
# \[s_{i+n-3} = \left( k_1+k_2-\sum_{j=3}^{i+2}k_j\right)^2, \quad \text{for}~~i=1,\dots,n-4\label{eq2}\]
# 
# \[s_{i+2n-7} = \left(\sum_{j=3}^{i+3} k_j\right)^2,\quad \text{for}~~i=1,\dots,n-4\label{eq3}\]
# 
# \[s =(k_1+k_2)^2\]
# 
# And those variables may depends on above variables but not linearly
# 
# \[(k_i\cdot k_j), \quad \text{for}~~i<j,~\text{and}~i,j\in\lbrace 4,\dots,n-1\rbrace\]
# 
# which can be used as \((n-4)(n-5)/2\) independent variables.
# 
# However, in total we have \(n(n-1)/2\) scalar products from
# \(k_1,k_2,\dots,k_n\).
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   We can take equations in \((\ref{eq1})\) by expansion
# 
#   \[s_i = k_1^2-2k_1\cdot\sum_{j=3}^{i+2}k_j +\left(\sum_{j=3}^{i+2}k_j\right)^2,\quad \text{for}~~i=1,\dots,n-3\]
# 
#   Then we can make subtraction \(s_i-s_{i-1}\) to obtain
# 
#   \[(k_1\cdot k_3) = \frac{1}{2}\left(k_1^2+k_3^2-s_1\right),\\(k_1\cdot k_4) = \frac{1}{2}\left(s_{1+2n-7}-k_3^2-s_2+s_1\right),\\(k_1\cdot k_{i+2}) = \frac{1}{2}\left(s_{i-1+2n-7}-s_{i-2+2n-7}-s_i+s_{i-1}\right),\quad \text{for}~~i=3,\dots,n-3\]
# 
#   And
# 
#   \[(k_1\cdot k_n) = k_1\cdot(k_1+k_2-k_3-\cdots-k_{n-1})\]
# 
#   Therefore we got \((n-2)\) scalar products
#   \((k_1\cdot k_3),~ (k_1\cdot p_4),~ \dots,~ (k_1\cdot k_n)\).
# \end{enumerate}
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   We can take equations in \((\ref{eq1})\) again by
# 
#   \[s_i = \left(k_2-\sum_{j=i+3}^n k_j\right)^2,\quad \text{for}~~i=1,\dots,n-3\]
# 
#   Then we can make subtraction \(s_i-s_{i+1}\) to obtain
# 
#   \[(k_2\cdot k_{i+3})=\frac{1}{2}(s_{i+n-3}-s_{i+1+n-3}-s_i+s_{i+1}),\quad \text{for}~~i=1,\dots,n-5\\(k_2\cdot k_{n-1})=\frac{1}{2}(s_{n-4+(n-3)}-k_n^2-s_i+s_{i+1}),\\(k_2\cdot k_n)=\frac{1}{2}(k_2^2+k_n^2-s_{n-3}),\]
# 
#   And according to \(s_{1+n-3}=(k_1+k_2-k_3)^2\) from \((\ref{eq2})\),
#   we can have
# 
#   \[(k_2\cdot k_3)=\frac{1}{2}(k_1^2+k_2^2+k_3^2+2k_1\cdot k_3-2k_1\cdot k_3-s_{1+n-3})\]
# 
#   Therefore, we got \(n-2\) scalar products
#   \((k_2\cdot k_3),~(k_2\cdot k_4),~\dots,~(k_2\cdot k_{n})\).
# \end{enumerate}
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   We can also take equations in \((\ref{eq3})\) by
# 
#   \[s_{i+2n-7} = \left(k_3+\sum_{j=4}^{i+3} k_j\right)^2,\quad \text{for}~~ i=1,\dots,n-4\]
# 
#   Then subtraction \(s_{i+2n-7}-s_{i-1+2n-7}\) can give us
# 
#   \[(k_3\cdot k_{i+3})=\frac{1}{2}\left[s_{i+2n-7}-s_{i-1+2n-7}-\left(\sum_{j=4}^{i+3}k_j\right)^2-\left(\sum_{j=4}^{i+2}k_j\right)^2\right],\quad \text{for}~~i=2,\dots,n-4\\(k_3\cdot k_4)=\frac{1}{2}(s_{1+(2n-7)}-k_3^2-k_4^2).\]
# 
#   And
# 
#   \[(k_3\cdot k_n) = k_3\cdot(k_1+k_2-k_3-\cdots-k_{n-1})\]
# 
#   Therefore, we got \(n-3\) scalar products
#   \((k_3\cdot k_4),~(k_3\cdot k_5),~\dots,~(k_3\cdot k_n)\).
# \item
#   We can take equations in \((\ref{eq2})\) by
# 
#   \[s_{i+(n-3)} = \left( k_n+\sum_{j=i+3}^{n-1}k_j\right)^2,\quad \text{for}~~i=1,\dots,n-4\]
# 
#   Then subtraction \(s_{i+(n-3)}-s_{i+1+(n-3)}\) can give us
# 
#   \[(k_n\cdot k_{i+3}) = \frac{1}{2}\left[s_{i+(n-3)}-s_{i+1+(n-3)}-\left(\sum_{j=i+3}^{n-1}k_j\right)^2+\left(\sum_{j=i+4}^{n-1}k_j\right)^2\right],\quad \text{for}~~i=1,\dots,n-5\\(k_n\cdot k_{n-1}) =\frac{1}{2}\left(s_{n-4+(n-3)}-k_n^2-k_{n-1}^2\right)\]
# 
#   Therefore we got \(n-4\) scalar products
#   \((k_n\cdot k_4),~(k_n\cdot k_5),~\dots,~(k_n\cdot k_{n-1})\).
# \end{enumerate}
# 
# \begin{enumerate}
# \def\labelenumi{\arabic{enumi}.}
# \item
#   Also we have \((k_1\cdot k_2) = \frac{1}{2} (s-k_1^2-k_2^2)\).
# \end{enumerate}
# 
# Finally in total we have
# \((n-2)+(n-2)+(n-3)+(n-4)+1+(n-4)(n-5)/2=n(n-1)/2\) solved scalar
# products, which is consistent with the total number.
# 
# \textbf{For the 1-(n-1) decay mode, we can simply change the sign of
# \(k_2\).}
#
##################################################################################################
"""
    generate_kin_relation( 
        n_inc::Int64, 
        n_out::Int64, 
        mom::Vector{Basic}, 
        mass2::Vector{Basic} 
    )::Dict{Basic,Basic}

Generate the kinematic relations, e.g. Mandelstam variables, according to the external fields.
"""
function generate_kin_relation( 
    n_inc::Int64, 
    n_out::Int64, 
    mom::Vector{Basic}, 
    mass2::Vector{Basic} 
)::Dict{Basic,Basic}
##################################################################################################

  @funs SP
  @vars shat

  nn = n_inc+n_out
  k2_sign = n_inc == 2 ? (+1) : (-1)
  half = Basic(1)/Basic(2)

  kin_relation = Dict{Basic,Basic}()

  ver_index_pre = 3*(nn-2)-4-1 # pre-occupied index for 3(n-2)-4 variable and one of them is shat
  ver_index = ver_index_pre + 1 
  for ii in 4:(nn-1)
    for jj in (ii+1):(nn-1)
      push!( kin_relation, make_SP(mom[ii],mom[jj]) => Basic("ver$(ver_index)") )
      ver_index += 1
    end # for jj
  end # for ii

  # on-shell conditions 
  for ii in 1:nn
    push!( kin_relation, make_SP(mom[ii],mom[ii]) => mass2[ii] )
  end # for ii


  # (k_1\cdot k_2) = \frac{1}{2} (s-k_1^2-k_2^2)  
  if n_out == 1
    push!( kin_relation, make_SP(mom[1],mom[2]) => half*( mass2[3] - mass2[1] - mass2[2] ) )
  else 
    push!( kin_relation, make_SP(mom[1],mom[2]) => k2_sign*half*( shat - mass2[1] - mass2[2] ) )
  end # if

  # (k_1\cdot k_3) = \frac{1}{2}\left(k_1^2+k_3^2-s_1\right),
  @vars ver1 # s_1
  if n_out == 1
    push!( kin_relation, make_SP(mom[1],mom[3]) => half*( mass2[1] + mass2[3] - mass2[2] ) )
  else 
    push!( kin_relation, make_SP(mom[1],mom[3]) => half*( mass2[1] + mass2[3] - ver1 ) )
  end # if

  @vars ver2 # s_2
  if nn >= 4
    if nn == 4
      # (k_1\cdot k_4) = \frac{1}{2}\left(shat-k_3^2-k_2^2+s_1\right),
      push!( kin_relation, make_SP(mom[1],mom[4]) => half*( shat - mass2[3] - mass2[2] + ver1 ) )
    else 
      # (k_1\cdot k_4) = \frac{1}{2}\left(s_{1+2n-7}-k_3^2-s_2+s_1\right),
      push!( kin_relation, make_SP(mom[1],mom[4]) => half*( Basic("ver$(1+2*nn-7)") - mass2[3] - ver2 + ver1 ) )
    end # if
  end # if

  # (k_1\cdot k_{i+2}) = \frac{1}{2}\left(s_{i-1+2n-7}-s_{i-2+2n-7}-s_i+s_{i-1}\right),\quad \text{for}~~i=3,\dots,n-3
  for ii in 3:(nn-3)
    push!( kin_relation, make_SP(mom[1],mom[ii+2]) => half*Basic("ver$(ii-1+2*nn-7) - ver$(ii-2+2*nn-7) - ver$(ii) + ver$(ii-1)") )
  end # for ii

  if nn > 4
    # (k_1\cdot k_n) = k_1\cdot(k_1+k_2-k_3-\cdots-k_{n-1})
    rhs = mass2[1] + k2_sign*make_SP(mom[1],mom[2]) # k_1\cdot k_1 + k_1\cdot k_2
    for ii = 3:(nn-1)
      rhs += (-1)*make_SP(mom[1],mom[ii])
    end # for ii
    rhs = subs( rhs, kin_relation )
    push!( kin_relation, make_SP(mom[1],mom[nn]) => rhs )
  end # if



  # (k_2\cdot k_{i+3})=\frac{1}{2}(s_{i+n-3}-s_{i+1+n-3}-s_i+s_{i+1}),\quad \text{for}~~i=1,\dots,n-5
  for ii in 1:(nn-5)
    push!( kin_relation, make_SP(mom[2],mom[ii+3]) => k2_sign*half*Basic("ver$(ii+nn-3) - ver$(ii+1+nn-3) - ver$(ii) + ver$(ii+1)") )
  end # for ii

  # (k_2\cdot k_n)=\frac{1}{2}(k_2^2+k_n^2-s_{n-3}),
  if n_out == 1
    push!( kin_relation, make_SP(mom[2],mom[3]) => half*(mass2[3]+mass2[2]-mass2[1]) )
  else 
    push!( kin_relation, make_SP(mom[2],mom[nn]) => k2_sign*half*(mass2[2]+mass2[nn]-Basic("ver$(nn-3)")) )
  end # if

  if nn >= 4
    if nn == 4
      # (k_2\cdot k_3) = k_1\cdot k_2 + k_2^2 - k_2\cdot k_4
      push!( kin_relation, make_SP(mom[2],mom[3]) => k2_sign*(expand∘subs)( k2_sign*make_SP(mom[1],mom[2]) + mass2[2] - k2_sign*make_SP(mom[2],mom[4]), kin_relation ) )
    else 
      # (k_2\cdot k_{n-1})=\frac{1}{2}(s_{n-4+(n-3)}-k_n^2-s_{n-4}+s_{n-3}),
      push!( kin_relation, make_SP(mom[2],mom[nn-1]) => k2_sign*half*Basic("ver$(nn-4+nn-3) - $(mass2[nn]) - ver$(nn-4) + ver$(nn-3)") )
    end # if
  end # if

  if nn == 4
    # (k_2\cdot k_3)=\frac{1}{2}(k_1^2+k_2^2+k_3^2+2k_1\cdot k_2-2k_1\cdot k_3-k_4^2)
    push!( kin_relation, make_SP(mom[2],mom[3]) => k2_sign*half*(expand∘subs)( mass2[1] + mass2[2] + mass2[3] + k2_sign*2*make_SP(mom[1],mom[2]) - 2*make_SP(mom[1],mom[3]) - mass2[4], kin_relation ) )
  elseif n_out == 1
    push!( kin_relation, make_SP(mom[2],mom[3]) => half*( mass2[3] + mass2[2] - mass2[1] ) )
  else 
    # (k_2\cdot k_3)=\frac{1}{2}(k_1^2+k_2^2+k_3^2+2k_1\cdot k_2-2k_1\cdot k_3-s_{1+n-3})
    push!( kin_relation, make_SP(mom[2],mom[3]) => k2_sign*half*(expand∘subs)( mass2[1] + mass2[2] + mass2[3] + k2_sign*2*make_SP(mom[1],mom[2]) - 2*make_SP(mom[1],mom[3]) - Basic("ver$(1+nn-3)"), kin_relation ) )
  end # if


  if nn >= 4
    if nn == 4
      # (k_3\cdot k_4)=\frac{1}{2}(shat-k_3^2-k_4^2).
      push!( kin_relation, make_SP(mom[3],mom[4]) => half*( shat - mass2[3] - mass2[4] ) )
    else 
      # (k_3\cdot k_4)=\frac{1}{2}(s_{1+(2n-7)}-k_3^2-k_4^2).
      push!( kin_relation, make_SP(mom[3],mom[4]) => half*( Basic("ver$(1+2*nn-7)") - mass2[3] - mass2[4] ) )
    end # if
  end # if


  # (k_3\cdot k_{i+3})=\frac{1}{2}\left[s_{i+2n-7}-s_{i-1+2n-7}-\left(\sum_{j=4}^{i+3}k_j\right)^2-\left(\sum_{j=4}^{i+2}k_j\right)^2\right],\quad \text{for}~~i=2,\dots,n-4
  for ii in 2:(nn-4)
    rhs = Basic("ver$(ii+2*nn-7) - ver$(ii-1+2*nn-7)") - mass2[ii+3]
    for jj in 4:(ii+2)
      rhs += (-2)*make_SP(mom[jj],mom[ii+3])
    end # for jj
    push!( kin_relation, make_SP(mom[3],mom[ii+3]) => half*(expand∘subs)( rhs, kin_relation ) )
  end # for ii

  # (k_3\cdot k_n) = k_3\cdot(k_1+k_2-k_3-\cdots-k_{n-1})
  rhs = make_SP(mom[1],mom[3])+k2_sign*make_SP(mom[2],mom[3])-mass2[3]
  for ii in 4:(nn-1)
    rhs += (-1)*make_SP(mom[3],mom[ii])
  end # for ii
  if n_out == 1
    push!( kin_relation, make_SP(mom[3],mom[3]) => mass2[3] )
  else
    push!( kin_relation, make_SP(mom[3],mom[nn]) => (expand∘subs)( rhs, kin_relation ) )
  end # if


  if nn >= 4
    if nn == 4 
      # (k_n\cdot k_{n-1}) =\frac{1}{2}\left(shat-k_n^2-k_{n-1}^2\right)
      push!( kin_relation, make_SP(mom[nn-1],mom[nn]) => half*( shat - mass2[nn] - mass2[nn-1] ) )
    else 
      # (k_n\cdot k_{n-1}) =\frac{1}{2}\left(s_{n-4+(n-3)}-k_n^2-k_{n-1}^2\right)
      push!( kin_relation, make_SP(mom[nn-1],mom[nn]) => half*( Basic("ver$(nn-4+nn-3)") - mass2[n] - mass2[nn-1] ) )
    end # if
  end # if

  # (k_n\cdot k_{i+3}) = \frac{1}{2}\left[s_{i+(n-3)}-s_{i+1+(n-3)}-\left(\sum_{j=i+3}^{n-1}k_j\right)^2+\left(\sum_{j=i+4}^{n-1}k_j\right)^2\right],\quad \text{for}~~i=1,\dots,n-5
  for ii in 1:(nn-5)
    rhs = Basic("ver$(ii+nn-3) - ver$(ii+1+nn-3)") - mass2[ii+3]
    for jj in (ii+4):(nn-1)
      rhs += (-2)*make_SP(mom[ii+3],mom[jj])
    end # for jj
    push!( kin_relation, make_SP(mom[ii+3],mom[nn]) => half*(expand∘subs)( rhs, kin_relation ) )
  end # for ii


  return kin_relation

end # function generate_kin_relation

















########################################################################################
"""
    generate_kin_relation( 
        n_inc::Int64, 
        n_out::Int64, 
        mom::Vector{Basic}, 
        mass2::Vector{Basic} 
    )::Dict{Basic,Basic}

"""
# Following the argument on the website:
#    https://9to5science.com/mandelstam-variables-for-2-to-3-particle-scattering
# 
# Generate the kinematic relations, e.g. Mandelstam variables, according to the external fields.
# 
# Consider all momentum are out-going, 
#   the Lorentz invariants for n-external-fields are defined as 
#     s^\prime_{i_1 \cdots i_\ell} = (p_{i_1}+\cdots+p_{i_\ell})^2.
# However, the expansion of the RHS shows that all of the s_{i_1 \cdots i_\ell} can be expressed 
#   as the linear combination of p_i\cdot p_j.
# Therefore, the Lorentz invariants for n-external-fields are 
#   s_{ij} = p_i\cdot p_j
# The counting number is C_n^2 = n(n-1)/2. 
# 
# Furthermore, we have the n-equation of nullity
#   \sum_{j=1;~j\ne i}^{n} s_{ij} = -m_i^2.
# So we need to remove n invariants, and the remaining counting number is
#   n(n-1)/2-n = n(n-3)/2.
# Specifically we remove the scalar products
#   p_1\cdot p_n, \dots, p_{n-1}\cdot p_n, and p_{n-2}\cdot p_{n-1}. 
# Finally the independent invariants are 
#   p_1\cdot p_2, \dots, p_1\cdot p_{n-1}, 
#   p_2\cdot p_3, \dots, p_2\cdot p_{n-1}, 
#     \cdots
#   p_{n-3}\cdot p_{n-2}, p_{n-3}\cdot p_{n-1}. 

function generate_kin_relation_v2( 
    n_inc::Int64, 
    n_out::Int64, 
    mom::Vector{Basic}, 
    mass2::Vector{Basic} 
)::Dict{Basic,Basic}
########################################################################################

  @funs SP
  @vars shat
  half = Basic(1)/Basic(2)

  nn = n_inc+n_out

  @assert n_inc in [1,2]
  mom_n = n_inc == 1 ?  mom[1] - sum(mom[2:nn-1]) : mom[1]+mom[2] - sum(mom[3:nn-1])

  # Assuming momenta are all outgoing, 
  k2_sign = n_inc == 2 ? (-1) : (+1) 
  # k1 is always incoming
  sign_list = Int64[ -1, k2_sign, fill(1,nn-2)... ]
  pp_list = sign_list .* mom
  #pp_list = [ -mom[1], k2_sign*mom[2], mom[3:nn]... ]

  kin_relation = Dict{Basic,Basic}()

  #------------
  if n_inc == 1 && n_out == 1
    push!( kin_relation, make_SP(mom[1],mom[2]) => mass2[1] )
  elseif n_inc == 2 && n_out == 1
    push!( kin_relation, make_SP(mom[1],mom[2]) => half*( mass2[3] - mass2[1] - mass2[2] ) )
  elseif n_inc == 1 && n_out == 2
    push!( kin_relation, make_SP(mom[1],mom[2]) => half*( mass2[1] + mass2[2] - mass2[3] ) )
  elseif n_inc == 2 && n_out >= 2
    push!( kin_relation, make_SP(mom[1],mom[2]) => half*( shat - mass2[1] - mass2[2] ) )
  end # if
  #------------


  ver_index = 1
  for ii in 1:(nn-3)
    for jj in (ii+1):(nn-1)
      if ii == 1 && jj == 2
        continue
      end # if 

      # v_ij = (p_i+p_j)^2, note here p_i may have minus sign compared to k_i
      push!( kin_relation, make_SP(mom[ii],mom[jj]) => sign_list[ii]*sign_list[jj]*half*( Basic("ver$(ver_index)") - mass2[ii] - mass2[jj] ) )
      ver_index += 1
    end # for jj
  end # for ii

  # on-shell conditions 
  for ii in 1:nn
    push!( kin_relation, make_SP(mom[ii],mom[ii]) => mass2[ii] )
  end # for ii

  # Then solve p_1\cdot p_n, \dots, p_{n-3}\cdot p_n via
  #   \sum_{j=1;~j\ne i}^{n} s_{ij} = -m_i^2 for i = 1,...,n-3
  for ii = 1:(nn-3)
    sum_pi_not_n = sum(pp_list[1:nn-1])-pp_list[ii]
    SP_expr = (expand∘subs)( -mass2[ii] - make_SP( pp_list[ii], sum_pi_not_n ), kin_relation )
    push!( kin_relation, make_SP(mom[ii],mom[nn]) => expand(sign_list[ii]*SP_expr) )
  end # for ii

  if nn >= 3
    # solve p_{n-2}\cdot p_{n-1} via 
    #   \sum_{i=1}^{n-1}\sum_{j=1}^{n-1} p_i\cdot p_j = m_n^2 
    SP_expr = subs( half*(mass2[nn] - make_SP( sum(pp_list[1:nn-1]), sum(pp_list[1:nn-1]) ) + 2*make_SP(pp_list[nn-2],pp_list[nn-1])), kin_relation )
    push!( kin_relation, make_SP( mom[nn-2], mom[nn-1] ) => expand( sign_list[nn-2]*sign_list[nn-1]*SP_expr) )
  end # if

  if nn >= 3
    # solve p_n\cdot p_{n-2} via
    #   \sum_{j=1;~j\ne n-2}^{n} s_{n-2,j} = -m_{n-2}^2 
    SP_expr = subs( -mass2[nn-2] - make_SP( pp_list[nn-2], sum(pp_list[1:nn-3]) ) - make_SP(pp_list[nn-2],pp_list[nn-1]), kin_relation )
    push!( kin_relation, make_SP( mom[nn], mom[nn-2] ) => expand( sign_list[nn]*sign_list[nn-2]*SP_expr ) )

    # solve p_n\cdot p_{n-1} via
    #   \sum_{j=1;~j\ne n-1}^{n} s_{n-1,j} = -m_{n-1}^2 
    SP_expr = subs( -mass2[nn-1] - make_SP( pp_list[nn-1], sum(pp_list[1:nn-2]) ), kin_relation )
    push!( kin_relation, make_SP( mom[nn], mom[nn-1] ) => expand( sign_list[nn]*sign_list[nn-1]*SP_expr ) )
  end # if

  return kin_relation

end # function generate_kin_relation_v2






########################################################################################
"""
    generate_kin_relation( graph_list::Vector{Graph} )::Dict{Basic,Basic}

Generate the kinematic relations for this processes including the internal non-loop propagators.
"""
function generate_kin_relation( graph_list::Vector{Graph} )::Dict{Basic,Basic}
########################################################################################

  graph0 = first( graph_list )

  n_inc = graph0.property[:n_inc]
  n_out = graph0.property[:n_out]
  nn = n_inc+n_out

  ext_edge_list = filter( e_ -> ( e_.property[:style]=="External" ), graph0.edge_list )
  @assert n_inc+n_out == length(ext_edge_list)

  sorted_ext_edge_list = sort( ext_edge_list, by=e_->e_.property[:mark] )
  mom = map( e_->e_.property[:momentum], sorted_ext_edge_list ) 
  mass2 = map( e_->e_.property[:particle].mass^2, sorted_ext_edge_list )
  # here mom and mass2 in fact are Vector{Basic}, but for later convenience we do not explicitly show it.

  @info "External Momenta" string(mom)

  kin_relation = generate_kin_relation_v2( n_inc, n_out, mom, mass2 )

  #--------------------------------------------------------
  # Next we need to define the denominators of propagators.
  #--------------------------------------------------------
  mom_n = Basic(0)
  for ii in 1:n_inc
    mom_n += mom[ii]
  end # for ii
  for ii in (n_inc+1):(nn-1)
    mom_n += (-1)*mom[ii]
  end # for ii

  @info "Momentum Conservation" "$(mom[nn]) = $(mom_n)"


  @funs Den
  den_set = Basic[]
  for g in graph_list
    int_edge_list = filter( e_->e_.property[:style] == "Internal", g.edge_list )
    for edge in int_edge_list
      den_mom = subs( edge.property[:momentum], mom[nn] => mom_n )
      den_mass = edge.property[:particle].mass
      den_width = edge.property[:particle].width
      union!( den_set, Den(expand(den_mom),den_mass,den_width) )
    end # for edge
  end # for g
  
  ext_mom_list = get_ext_momenta(den_set)
  sort!(den_set, by=den -> begin
      tmp_mom, tmp_mass, tmp_width = get_args(den)
      if SymEngine.coeff(tmp_mom, first(ext_mom_list)) < 0
        tmp_mom = -tmp_mom
      end # if
      tmp_mom = expand(tmp_mom)
      tmp_mass = expand(tmp_mass)
      tmp_width = expand(tmp_width)
      tmp_ext_coeff_mat = (vec∘coefficient_matrix)( [tmp_mom], ext_mom_list )
      tmp_k_str = (join∘map)( coeff->coeff==-1 ? "2" : string(coeff), tmp_ext_coeff_mat)
      (
        parse(BigInt, reverse(tmp_k_str), base=3),
        string(tmp_mass),
        string(tmp_width)
      )
    end
  ) # end sort!

  symbol_list = (free_symbols∘collect∘values)(kin_relation)
  ver_list = filter( x-> (length∘string)(x)>3&&string(x)[1:3]=="ver", symbol_list )
  if nn >= 3
    @assert length(ver_list) in [ nn*(nn-3)/Basic(2), nn*(nn-3)/Basic(2)-one(Basic) ]
  end # if

  # we know at most n(n-3)/2 verI's have been occupied
  ver_index_pre = length(ver_list) 
  ver_index = ver_index_pre + 1
  for one_den in den_set
    @assert Basic("ver$(ver_index)") ∉ symbol_list 

    arg_list = get_args(one_den)
    den_mom = arg_list[1]
    den_mass = arg_list[2]
    den_width = arg_list[3]
    push!( kin_relation, Den(expand(den_mom),den_mass,den_width) => Basic("ver$(ver_index)") )
    push!( kin_relation, Den(expand(-den_mom),den_mass,den_width) => Basic("ver$(ver_index)") )

    ver_index += 1 
  end # for one_den

  for g in graph_list
    int_edge_list = filter( e_->e_.property[:style] == "Internal", g.edge_list )
    for edge in int_edge_list
      den_mom = edge.property[:momentum]
      den_mass = edge.property[:particle].mass
      den_width = edge.property[:particle].width

      subs_den_mom = (expand∘subs)( den_mom, mom[nn] => mom_n )
      if (!iszero∘expand)( den_mom - subs_den_mom )
        push!( kin_relation, Den(expand(den_mom),den_mass,den_width) => subs(Den(subs_den_mom,den_mass,den_width),kin_relation) )
        push!( kin_relation, Den(expand(-den_mom),den_mass,den_width) => subs(Den(subs_den_mom,den_mass,den_width),kin_relation) )
      end # if

    end # for edge
  end # for g

  return kin_relation

end # function generate_kin_relation

#################################################################################
"""
    generate_ext_mom_list( graph_list::Vector{Graph} )::Vector{Basic}

Generate the list of external momenta according to this process.
"""
function generate_ext_mom_list( graph_list::Vector{Graph} )::Vector{Basic}
#################################################################################

  graph0 = first( graph_list )

  n_inc = graph0.property[:n_inc]
  n_out = graph0.property[:n_out]

  ext_edge_list = filter( e_ -> ( e_.property[:style]=="External" ), graph0.edge_list )
  @assert n_inc+n_out == length(ext_edge_list)

  sorted_ext_edge_list = sort( ext_edge_list, by=e_->e_.property[:mark] )
  ext_mom_list = map( e_->e_.property[:momentum], sorted_ext_edge_list ) 

  return ext_mom_list

end # function generate_ext_mom_list






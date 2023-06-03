# [Functions in Kin.jl](@id FunctionsKin)

```@docs
FeAmGen.generate_gauge_choice
FeAmGen.generate_kin_relation
FeAmGen.generate_ext_mom_list
```

### Kinematic variables for (``n-2``)-body final state.

Considering ``2\to (n-2)`` process, we set up momenta
``k_1,k_2,k_3,\dots,k_n``.

### According to RAJENDRA KUMAR, Phys.Rev. D2 (1970) 1902-1914

For (``n-2``)-body final state, there Eq. (2.6) shows the ``(3(n-2)-4)``
independent variables are

```math
s_{i} = \left( k_1 -\sum_{j=3}^{i+2} k_j \right)^2, \quad \text{for}~~i=1,\dots,n-3
\qquad (1)
```

```math
s_{i+n-3} = \left( k_1+k_2-\sum_{j=3}^{i+2}k_j\right)^2, \quad \text{for}~~i=1,\dots,n-4
\qauad (2)
```

```math
s_{i+2n-7} = \left(\sum_{j=3}^{i+3} k_j\right)^2,\quad \text{for}~~i=1,\dots,n-4 
\qauad (3)
```

```math
s =(k_1+k_2)^2 
```

And those variables may depends on above variables but not linearly

```math
(k_i\cdot k_j), \quad \text{for}~~i<j,~\text{and}~i,j\in\lbrace 4,\dots,n-1\rbrace
```

which can be used as ``(n-4)(n-5)/2`` independent variables.

However, in total we have ``n(n-1)/2`` scalar products from
``k_1,k_2,\dots,k_n``.

+ We can take equations in (1) by expansion

  ```math
  s_i = k_1^2-2k_1\cdot\sum_{j=3}^{i+2}k_j +\left(\sum_{j=3}^{i+2}k_j\right)^2,\quad \text{for}~~i=1,\dots,n-3
  ```

  Then we can make subtraction ``s_i-s_{i-1}`` to obtain

  ```math
  (k_1\cdot k_3) = \frac{1}{2}\left(k_1^2+k_3^2-s_1\right),\\(k_1\cdot k_4) = \frac{1}{2}\left(s_{1+2n-7}-k_3^2-s_2+s_1\right),\\(k_1\cdot k_{i+2}) = \frac{1}{2}\left(s_{    i-1+2n-7}-s_{i-2+2n-7}-s_i+s_{i-1}\right),\quad \text{for}~~i=3,\dots,n-3
  ```

  And

  ``(k_1\cdot k_n) = k_1\cdot(k_1+k_2-k_3-\cdots-k_{n-1})``

  Therefore we got ``(n-2)`` scalar products
  ``(k_1\cdot k_3),~ (k_1\cdot p_4),~ \dots,~ (k_1\cdot k_n)``.


+ We can also take equations in (3) by

  ```math
  s_{i+2n-7} = \left(k_3+\sum_{j=4}^{i+3} k_j\right)^2,\quad \text{for}~~ i=1,\dots,n-4
  ```

  Then subtraction ``s_{i+2n-7}-s_{i-1+2n-7}`` can give us

  ```math
   (k_3\cdot k_{i+3})=\frac{1}{2}\left[s_{i+2n-7}-s_{i-1+2n-7}-\left(\sum_{j=4}^{i+3}k_j\right)^2-\left(\sum_{j=4}^{i+2}k_j\right)^2\right],\quad \text{for}~~i=2,\dots,n-4    \\(k_3\cdot k_4)=\frac{1}{2}(s_{1+(2n-7)}-k_3^2-k_4^2).
  ```

  And

  ```math
  (k_3\cdot k_n) = k_3\cdot(k_1+k_2-k_3-\cdots-k_{n-1})
  ```

  Therefore, we got ``n-3`` scalar products
  ```math
  (k_3\cdot k_4),~(k_3\cdot k_5),~\dots,~(k_3\cdot k_n).
  ```

+  We can take equations in (2) by

  ```math
  s_{i+(n-3)} = \left( k_n+\sum_{j=i+3}^{n-1}k_j\right)^2,\quad \text{for}~~i=1,\dots,n-4
  ```

  Then subtraction ``s_{i+(n-3)}-s_{i+1+(n-3)}`` can give us

  ```math
  (k_n\cdot k_{i+3}) = \frac{1}{2}\left[s_{i+(n-3)}-s_{i+1+(n-3)}-\left(\sum_{j=i+3}^{n-1}k_j\right)^2+\left(\sum_{j=i+4}^{n-1}k_j\right)^2\right],\quad \text{for}~~i=1,\dots,n-5\\(k_n\cdot k_{n-1}) =\frac{1}{2}\left(s_{n-4+(n-3)}-k_n^2-k_{n-1}^2\right)
  ```

  Therefore we got ``n-4`` scalar products
  ```math
  (k_n\cdot k_4),~(k_n\cdot k_5),~\dots,~(k_n\cdot k_{n-1}).
  ```

+  Also we have ``(k_1\cdot k_2) = \frac{1}{2} (s-k_1^2-k_2^2)``.



Finally in total we have
``(n-2)+(n-2)+(n-3)+(n-4)+1+(n-4)(n-5)/2=n(n-1)/2`` solved scalar
products, which is consistent with the total number.

### For the 1-(n-1) decay mode, we can simply change the sign of ``k_2``.





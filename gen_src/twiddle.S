#define       Z              %rdi
#define       W              %rsi
#define       k2             %rdx
#define       N              %rcx

              .global        twiddle_simd4
              .global        twiddle_simd2
              .global        twiddle_scalar

twiddle_simd4:
              xor            %rax, %rax
loop1:        mov            %rax, %r9
              mov            %rax, %r8
              shl            %r9
              imul           k2, %r8
              vbroadcastsd   (W, %r8, 8), %ymm0
              vbroadcastsd   8(W, %r8, 8), %ymm1
              vmulpd         (Z, %r9, 8), %ymm0, %ymm2
              vmulpd         32(Z, %r9, 8), %ymm0, %ymm3
              vfmsub231pd    32(Z, %r9, 8), %ymm1, %ymm2
              vfmadd231pd    (Z, %r9, 8), %ymm1, %ymm3
              vmovapd        %ymm2, (Z, %r9, 8)
              vmovapd        %ymm3, 32(Z, %r9, 8)
              cmp            %rax, N
              jne            loop1
              ret

twiddle_simd2:
              xor            %rax, %rax
loop2:        mov            %rax, %r8
              shl            %r9
              imul           k2, %r8
              vpbroadcastq   (W, %r8, 8), %xmm0
              vpbroadcastq   8(W, %r8, 8), %xmm1
              vmulpd         (Z, %r9, 8), %xmm0, %xmm2
              vmulpd         16(Z, %r9, 8), %xmm0, %xmm3
              vfmsub231pd    16(Z, %r9, 8), %xmm1, %xmm2
              vfmadd231pd    (Z, %r9, 8), %xmm1, %xmm3
              vmovapd        %xmm2, (Z, %r9, 8)
              vmovapd        %xmm3, 32(Z, %r9, 8)
              cmp            %rax, N
              jne            loop2
              ret

twiddle_scalar:
              xor            %rax, %rax
loop3:        mov            %rax, %r8
              shl            %r9
              imul           k2, %r8
              vmovq          (W, %r8, 8), %xmm0
              vmovq          8(W, %r8, 8), %xmm1
              vmulpd         (Z, %r9, 8), %xmm0, %xmm2
              vmulpd         8(Z, %r9, 8), %xmm0, %xmm3
              vfmsub231pd    8(Z, %r9, 8), %xmm1, %xmm2
              vfmadd231pd    (Z, %r9, 8), %xmm1, %xmm3
              vmovapd        %xmm2, (Z, %r9, 8)
              vmovapd        %xmm3, 32(Z, %r9, 8)
              cmp            %rax, N
              jne            loop3
              ret


#define       Z              %rdi
#define       W              %rsi
#define       k2             %rdx
#define       N              %rcx

              .global        sfft_twiddle_w4
              .global        sfft_twiddle_w2
              .global        sfft_twiddle_w1

              .text

sfft_twiddle_w4:
              xor            %rax, %rax
loop1:        mov            %rax, %r9
              mov            %rax, %r8
              shl            $3, %r9
              shl            %r8
              imul           k2, %r8
              vbroadcastsd   (W, %r8, 8), %ymm0
              vbroadcastsd   8(W, %r8, 8), %ymm1
              vmulpd         32(Z, %r9, 8), %ymm1, %ymm2
              vmulpd         (Z, %r9, 8), %ymm1, %ymm3
              vfmsub231pd    (Z, %r9, 8), %ymm0, %ymm2
              vfmadd231pd    32(Z, %r9, 8), %ymm0, %ymm3
              vmovapd        %ymm2, (Z, %r9, 8)
              vmovapd        %ymm3, 32(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop1
              ret

sfft_twiddle_w2:
              xor            %rax, %rax
loop2:        mov            %rax, %r9
              mov            %rax, %r8
              shl            $2, %r9
              shl            %r8
              imul           k2, %r8
              vpbroadcastq   (W, %r8, 8), %xmm0
              vpbroadcastq   8(W, %r8, 8), %xmm1
              vmulpd         16(Z, %r9, 8), %xmm1, %xmm2
              vmulpd         (Z, %r9, 8), %xmm1, %xmm3
              vfmsub231pd    (Z, %r9, 8), %xmm0, %xmm2
              vfmadd231pd    16(Z, %r9, 8), %xmm0, %xmm3
              vmovapd        %xmm2, (Z, %r9, 8)
              vmovapd        %xmm3, 16(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop2
              ret

sfft_twiddle_w1:
              xor            %rax, %rax
loop3:        mov            %rax, %r9
              mov            %rax, %r8
              shl            %r9
              shl            %r8
              imul           k2, %r8
              vmovq          (W, %r8, 8), %xmm0
              vmovq          8(W, %r8, 8), %xmm1
              vmulsd         8(Z, %r9, 8), %xmm1, %xmm2
              vmulsd         (Z, %r9, 8), %xmm1, %xmm3
              vfmsub231sd    (Z, %r9, 8), %xmm0, %xmm2
              vfmadd231sd    8(Z, %r9, 8), %xmm0, %xmm3
              vmovq          %xmm2, (Z, %r9, 8)
              vmovq          %xmm3, 8(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop3
              ret


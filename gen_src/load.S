#define       X              %rdi
#define       Y              %rsi
#define       Z              %rdx
#define       stride         %rcx
#define       N              %r8

              .global        sfft_load_complex_w4
              .global        sfft_load_complex_w2
              .global        sfft_load_complex_w1

              .text

sfft_load_complex_w4:
              xor            %rax, %rax
loop1:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            $3, %r9
              vmovupd        (X, %r10, 8), %ymm0
              vmovupd        (Y, %r10, 8), %ymm1
              vmovupd        %ymm0, (Z, %r9, 8)
              vmovupd        %ymm1, 32(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop1
              ret

sfft_load_complex_w2:
              xor            %rax, %rax
loop2:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            $2, %r9
              vmovupd        (X, %r10, 8), %xmm0
              vmovupd        (Y, %r10, 8), %xmm1
              vmovupd        %xmm0, (Z, %r9, 8)
              vmovupd        %xmm1, 16(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop2
              ret

sfft_load_complex_w1:
              xor            %rax, %rax
loop3:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            %r9
              vmovq          (X, %r10, 8), %xmm0
              vmovq          (Y, %r10, 8), %xmm1
              vmovq          %xmm0, (Z, %r9, 8)
              vmovq          %xmm1, 8(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop3
              ret

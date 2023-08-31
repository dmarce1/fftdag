#define       X              %rdi
#define       Y              %rsi
#define       Z              %rdx
#define       stride         %rcx
#define       N              %r8

              .global        load_complex_simd4
              .global        load_complex_simd2
              .global        load_complex_scalar

load_complex_simd4:
              xor            %rax, %rax
loop1:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            %r9
              vmovapd        (X, %r10, 8), %ymm0
              vmovapd        (Y, %r10, 8), %ymm1
              vmovupd        %ymm0, (Z, %r9, 8)
              vmovupd        %ymm1, 32(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop1
              ret

load_complex_simd2:
              xor            %rax, %rax
loop2:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            %r9
              vmovapd        (X, %r10, 8), %xmm0
              vmovapd        (Y, %r10, 8), %xmm1
              vmovupd        %xmm0, (Z, %r9, 8)
              vmovupd        %xmm1, 16(Z, %r9, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop2
              ret

load_complex_scalar:
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

#define       X              %rdi
#define       Y              %rsi
#define       Z              %rdx
#define       stride         %rcx
#define       N              %r8

              .global        store_complex_simd4
              .global        store_complex_simd2
              .global        store_complex_scalar

store_complex_simd4:
              xor            %rax, %rax
loop1:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            %r9
              vmovupd        (Z, %r9, 8), %ymm0
              vmovupd        32(Z, %r9, 8), %ymm1
              vmovapd        %ymm0, (X, %r10, 8)
              vmovapd        %ymm1, (Y, %r10, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop1
              ret

store_complex_simd2:
              xor            %rax, %rax
loop2:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            %r9
              vmovupd        (Z, %r9, 8), %xmm0
              vmovupd        16(Z, %r9, 8), %xmm1
              vmovapd        %xmm0, (X, %r10, 8)
              vmovapd        %xmm1, (Y, %r10, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop2
              ret

store_complex_scalar:
              xor            %rax, %rax
loop3:        mov            stride, %r10
              mov            %rax, %r9
              imul           %rax, %r10
              shl            %r9
              vmovq          (Z, %r9, 8), %xmm0
              vmovq          8(Z, %r9, 8), %xmm1
              vmovq          %xmm0, (X, %r10, 8)
              vmovq          %xmm1, (Y, %r10, 8)
              inc            %rax
              cmp            %rax, N
              jne            loop3
              ret

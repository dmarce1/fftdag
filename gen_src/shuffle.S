              .global        sfft_perf_shuf_w4
              .global        sfft_perf_shuf_w2
              .global        sfft_perf_shuf_w1
              .global        sfft_inv_perf_shuf_w4
              .global        sfft_inv_perf_shuf_w2
              .global        sfft_inv_perf_shuf_w1

              .text

sfft_perf_shuf_w4:
              vmovupd        (%rdi), %ymm0
              vperm2f128     $0x20, (%rsi), %ymm0, %ymm2
              vperm2f128     $0x31, (%rsi), %ymm0, %ymm3
              vpermpd        $216, %ymm2, %ymm0
              vpermpd        $216, %ymm3, %ymm1
              vmovupd        %ymm0, (%rdi)
              vmovupd        %ymm1, (%rsi)
              ret

sfft_perf_shuf_w2:
              vmovq          8(%rdi), %xmm0
              vmovq          (%rsi), %xmm1
              vmovq          %xmm1, 8(%rdi)
              vmovq          %xmm0, (%rsi)
              ret

sfft_perf_shuf_w1:
              ret


sfft_inv_perf_shuf_w4:
              vpermpd        $216, (%rdi), %ymm0
              vpermpd        $216, (%rsi), %ymm1
              vperm2f128     $0x20, %ymm1, %ymm0, %ymm2
              vperm2f128     $0x31, %ymm1, %ymm0, %ymm3
              vmovupd        %ymm2, (%rdi)
              vmovupd        %ymm3, (%rsi)
              ret


sfft_inv_perf_shuf_w2:
              vmovq          8(%rdi), %xmm0
              vmovq          (%rsi), %xmm1
              vmovq          %xmm1, 8(%rdi)
              vmovq          %xmm0, (%rsi)
              ret

sfft_inv_perf_shuf_w1:
              ret


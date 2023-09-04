
              .global        sfft_complex_w4
              .global        sfft_complex_w2
              .global        sfft_complex_w1
              .global        sfft_complex_dit_w4
              .global        sfft_complex_dit_w2
              .global        sfft_complex_dit_w1
              .global        sfft_complex_dif_w4
              .global        sfft_complex_dif_w2
              .global        sfft_complex_dif_w1
              .text

sfft_complex_w1:
              lea            (functions1), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_w2:
              lea            (functions2), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_w4:
              lea            (functions4), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_dit_w1:
              lea            (functions5), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_dit_w2:
              lea            (functions6), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_dit_w4:
              lea            (functions7), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_dif_w1:
              lea            (functions8), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_dif_w2:
              lea            (functions9), %rax
              call           *(%rax, %rcx, 8)
		      ret

sfft_complex_dif_w4:
              lea            (functions10), %rax
              call           *(%rax, %rcx, 8)
		      ret

              .align         8
functions1:   .quad          0
              .quad          1
              .quad          sfft_complex_w1_2
              .quad          sfft_complex_w1_3
              .quad          sfft_complex_w1_4
              .quad          sfft_complex_w1_5
              .quad          sfft_complex_w1_6
              .quad          sfft_complex_w1_7
              .quad          sfft_complex_w1_8
              .quad          sfft_complex_w1_9
              .quad          sfft_complex_w1_10
              .quad          sfft_complex_w1_11
              .quad          sfft_complex_w1_12
              .quad          sfft_complex_w1_13
              .quad          sfft_complex_w1_14
              .quad          sfft_complex_w1_15
              .quad          sfft_complex_w1_16
              .quad          sfft_complex_w1_17
              .quad          sfft_complex_w1_18
              .quad          sfft_complex_w1_19
              .quad          sfft_complex_w1_20
              .quad          sfft_complex_w1_21
              .quad          sfft_complex_w1_22
              .quad          sfft_complex_w1_23
              .quad          sfft_complex_w1_24
              .quad          sfft_complex_w1_25
              .quad          sfft_complex_w1_26
              .quad          sfft_complex_w1_27
              .quad          sfft_complex_w1_28
              .quad          sfft_complex_w1_29
              .quad          sfft_complex_w1_30
functions2:   .quad          sfft_complex_w1_31
              .quad          sfft_complex_w1_32
              .quad          sfft_complex_w2_2
              .quad          sfft_complex_w2_3
              .quad          sfft_complex_w2_4
              .quad          sfft_complex_w2_5
              .quad          sfft_complex_w2_6
              .quad          sfft_complex_w2_7
              .quad          sfft_complex_w2_8
              .quad          sfft_complex_w2_9
              .quad          sfft_complex_w2_10
              .quad          sfft_complex_w2_11
              .quad          sfft_complex_w2_12
              .quad          sfft_complex_w2_13
              .quad          sfft_complex_w2_14
              .quad          sfft_complex_w2_15
              .quad          sfft_complex_w2_16
              .quad          sfft_complex_w2_17
              .quad          sfft_complex_w2_18
              .quad          sfft_complex_w2_19
              .quad          sfft_complex_w2_20
              .quad          sfft_complex_w2_21
              .quad          sfft_complex_w2_22
              .quad          sfft_complex_w2_23
              .quad          sfft_complex_w2_24
              .quad          sfft_complex_w2_25
              .quad          sfft_complex_w2_26
              .quad          sfft_complex_w2_27
              .quad          sfft_complex_w2_28
              .quad          sfft_complex_w2_29
              .quad          sfft_complex_w2_30
functions4:   .quad          sfft_complex_w2_31
              .quad          sfft_complex_w2_32
              .quad          sfft_complex_w4_2
              .quad          sfft_complex_w4_3
              .quad          sfft_complex_w4_4
              .quad          sfft_complex_w4_5
              .quad          sfft_complex_w4_6
              .quad          sfft_complex_w4_7
              .quad          sfft_complex_w4_8
              .quad          sfft_complex_w4_9
              .quad          sfft_complex_w4_10
              .quad          sfft_complex_w4_11
              .quad          sfft_complex_w4_12
              .quad          sfft_complex_w4_13
              .quad          sfft_complex_w4_14
              .quad          sfft_complex_w4_15
              .quad          sfft_complex_w4_16
              .quad          sfft_complex_w4_17
              .quad          sfft_complex_w4_18
              .quad          sfft_complex_w4_19
              .quad          sfft_complex_w4_20
              .quad          sfft_complex_w4_21
              .quad          sfft_complex_w4_22
              .quad          sfft_complex_w4_23
              .quad          sfft_complex_w4_24
              .quad          sfft_complex_w4_25
              .quad          sfft_complex_w4_26
              .quad          sfft_complex_w4_27
              .quad          sfft_complex_w4_28
              .quad          sfft_complex_w4_29
              .quad          sfft_complex_w4_30
              .quad          sfft_complex_w4_31
              .quad          sfft_complex_w4_32
functions5:   .quad          0
              .quad          1
              .quad          sfft_complex_dit_w1_2
              .quad          sfft_complex_dit_w1_3
              .quad          sfft_complex_dit_w1_4
              .quad          sfft_complex_dit_w1_5
              .quad          sfft_complex_dit_w1_6
              .quad          sfft_complex_dit_w1_7
              .quad          sfft_complex_dit_w1_8
              .quad          sfft_complex_dit_w1_9
              .quad          sfft_complex_dit_w1_10
              .quad          sfft_complex_dit_w1_11
              .quad          sfft_complex_dit_w1_12
              .quad          sfft_complex_dit_w1_13
              .quad          sfft_complex_dit_w1_14
              .quad          sfft_complex_dit_w1_15
              .quad          sfft_complex_dit_w1_16
              .quad          sfft_complex_dit_w1_17
              .quad          sfft_complex_dit_w1_18
              .quad          sfft_complex_dit_w1_19
              .quad          sfft_complex_dit_w1_20
              .quad          sfft_complex_dit_w1_21
              .quad          sfft_complex_dit_w1_22
              .quad          sfft_complex_dit_w1_23
              .quad          sfft_complex_dit_w1_24
              .quad          sfft_complex_dit_w1_25
              .quad          sfft_complex_dit_w1_26
              .quad          sfft_complex_dit_w1_27
              .quad          sfft_complex_dit_w1_28
              .quad          sfft_complex_dit_w1_29
              .quad          sfft_complex_dit_w1_30
functions6:   .quad          sfft_complex_dit_w1_31
              .quad          sfft_complex_dit_w1_32
              .quad          sfft_complex_dit_w2_2
              .quad          sfft_complex_dit_w2_3
              .quad          sfft_complex_dit_w2_4
              .quad          sfft_complex_dit_w2_5
              .quad          sfft_complex_dit_w2_6
              .quad          sfft_complex_dit_w2_7
              .quad          sfft_complex_dit_w2_8
              .quad          sfft_complex_dit_w2_9
              .quad          sfft_complex_dit_w2_10
              .quad          sfft_complex_dit_w2_11
              .quad          sfft_complex_dit_w2_12
              .quad          sfft_complex_dit_w2_13
              .quad          sfft_complex_dit_w2_14
              .quad          sfft_complex_dit_w2_15
              .quad          sfft_complex_dit_w2_16
              .quad          sfft_complex_dit_w2_17
              .quad          sfft_complex_dit_w2_18
              .quad          sfft_complex_dit_w2_19
              .quad          sfft_complex_dit_w2_20
              .quad          sfft_complex_dit_w2_21
              .quad          sfft_complex_dit_w2_22
              .quad          sfft_complex_dit_w2_23
              .quad          sfft_complex_dit_w2_24
              .quad          sfft_complex_dit_w2_25
              .quad          sfft_complex_dit_w2_26
              .quad          sfft_complex_dit_w2_27
              .quad          sfft_complex_dit_w2_28
              .quad          sfft_complex_dit_w2_29
              .quad          sfft_complex_dit_w2_30
functions7:   .quad          sfft_complex_dit_w2_31
              .quad          sfft_complex_dit_w2_32
              .quad          sfft_complex_dit_w4_2
              .quad          sfft_complex_dit_w4_3
              .quad          sfft_complex_dit_w4_4
              .quad          sfft_complex_dit_w4_5
              .quad          sfft_complex_dit_w4_6
              .quad          sfft_complex_dit_w4_7
              .quad          sfft_complex_dit_w4_8
              .quad          sfft_complex_dit_w4_9
              .quad          sfft_complex_dit_w4_10
              .quad          sfft_complex_dit_w4_11
              .quad          sfft_complex_dit_w4_12
              .quad          sfft_complex_dit_w4_13
              .quad          sfft_complex_dit_w4_14
              .quad          sfft_complex_dit_w4_15
              .quad          sfft_complex_dit_w4_16
              .quad          sfft_complex_dit_w4_17
              .quad          sfft_complex_dit_w4_18
              .quad          sfft_complex_dit_w4_19
              .quad          sfft_complex_dit_w4_20
              .quad          sfft_complex_dit_w4_21
              .quad          sfft_complex_dit_w4_22
              .quad          sfft_complex_dit_w4_23
              .quad          sfft_complex_dit_w4_24
              .quad          sfft_complex_dit_w4_25
              .quad          sfft_complex_dit_w4_26
              .quad          sfft_complex_dit_w4_27
              .quad          sfft_complex_dit_w4_28
              .quad          sfft_complex_dit_w4_29
              .quad          sfft_complex_dit_w4_30
              .quad          sfft_complex_dit_w4_31
              .quad          sfft_complex_dit_w4_32
functions8:   .quad          0
              .quad          1
              .quad          sfft_complex_dif_w1_2
              .quad          sfft_complex_dif_w1_3
              .quad          sfft_complex_dif_w1_4
              .quad          sfft_complex_dif_w1_5
              .quad          sfft_complex_dif_w1_6
              .quad          sfft_complex_dif_w1_7
              .quad          sfft_complex_dif_w1_8
              .quad          sfft_complex_dif_w1_9
              .quad          sfft_complex_dif_w1_10
              .quad          sfft_complex_dif_w1_11
              .quad          sfft_complex_dif_w1_12
              .quad          sfft_complex_dif_w1_13
              .quad          sfft_complex_dif_w1_14
              .quad          sfft_complex_dif_w1_15
              .quad          sfft_complex_dif_w1_16
              .quad          sfft_complex_dif_w1_17
              .quad          sfft_complex_dif_w1_18
              .quad          sfft_complex_dif_w1_19
              .quad          sfft_complex_dif_w1_20
              .quad          sfft_complex_dif_w1_21
              .quad          sfft_complex_dif_w1_22
              .quad          sfft_complex_dif_w1_23
              .quad          sfft_complex_dif_w1_24
              .quad          sfft_complex_dif_w1_25
              .quad          sfft_complex_dif_w1_26
              .quad          sfft_complex_dif_w1_27
              .quad          sfft_complex_dif_w1_28
              .quad          sfft_complex_dif_w1_29
              .quad          sfft_complex_dif_w1_30
functions9:   .quad          sfft_complex_dif_w1_31
              .quad          sfft_complex_dif_w1_32
              .quad          sfft_complex_dif_w2_2
              .quad          sfft_complex_dif_w2_3
              .quad          sfft_complex_dif_w2_4
              .quad          sfft_complex_dif_w2_5
              .quad          sfft_complex_dif_w2_6
              .quad          sfft_complex_dif_w2_7
              .quad          sfft_complex_dif_w2_8
              .quad          sfft_complex_dif_w2_9
              .quad          sfft_complex_dif_w2_10
              .quad          sfft_complex_dif_w2_11
              .quad          sfft_complex_dif_w2_12
              .quad          sfft_complex_dif_w2_13
              .quad          sfft_complex_dif_w2_14
              .quad          sfft_complex_dif_w2_15
              .quad          sfft_complex_dif_w2_16
              .quad          sfft_complex_dif_w2_17
              .quad          sfft_complex_dif_w2_18
              .quad          sfft_complex_dif_w2_19
              .quad          sfft_complex_dif_w2_20
              .quad          sfft_complex_dif_w2_21
              .quad          sfft_complex_dif_w2_22
              .quad          sfft_complex_dif_w2_23
              .quad          sfft_complex_dif_w2_24
              .quad          sfft_complex_dif_w2_25
              .quad          sfft_complex_dif_w2_26
              .quad          sfft_complex_dif_w2_27
              .quad          sfft_complex_dif_w2_28
              .quad          sfft_complex_dif_w2_29
              .quad          sfft_complex_dif_w2_30
functions10:  .quad          sfft_complex_dif_w2_31
              .quad          sfft_complex_dif_w2_32
              .quad          sfft_complex_dif_w4_2
              .quad          sfft_complex_dif_w4_3
              .quad          sfft_complex_dif_w4_4
              .quad          sfft_complex_dif_w4_5
              .quad          sfft_complex_dif_w4_6
              .quad          sfft_complex_dif_w4_7
              .quad          sfft_complex_dif_w4_8
              .quad          sfft_complex_dif_w4_9
              .quad          sfft_complex_dif_w4_10
              .quad          sfft_complex_dif_w4_11
              .quad          sfft_complex_dif_w4_12
              .quad          sfft_complex_dif_w4_13
              .quad          sfft_complex_dif_w4_14
              .quad          sfft_complex_dif_w4_15
              .quad          sfft_complex_dif_w4_16
              .quad          sfft_complex_dif_w4_17
              .quad          sfft_complex_dif_w4_18
              .quad          sfft_complex_dif_w4_19
              .quad          sfft_complex_dif_w4_20
              .quad          sfft_complex_dif_w4_21
              .quad          sfft_complex_dif_w4_22
              .quad          sfft_complex_dif_w4_23
              .quad          sfft_complex_dif_w4_24
              .quad          sfft_complex_dif_w4_25
              .quad          sfft_complex_dif_w4_26
              .quad          sfft_complex_dif_w4_27
              .quad          sfft_complex_dif_w4_28
              .quad          sfft_complex_dif_w4_29
              .quad          sfft_complex_dif_w4_30
              .quad          sfft_complex_dif_w4_31
              .quad          sfft_complex_dif_w4_32




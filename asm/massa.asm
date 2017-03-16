; ----------------------------------------------------------------------------------------
; Massa di una stella. Linux 64 bit
; 
;Assembla ed esegue:
;     nasm -felf64 massa.asm && gcc massa.o && ./a.out
; ----------------------------------------------------------------------------------------
        global  main
        extern  printf


        section .data
        
        n_range:
         TIMES 8 dq 0
        h:
         dq 0.1
        h_half:
         dq 0.05
        mesh:
         TIMES 100 dq 0
        theta:
         TIMES 100 dq 0
        phi:
         TIMES 100 dq 0
        half:
         dq 0.5
        zero:
         dq 0.0
        uno:
         dq 1.0
        meno_due:
         dq -2.0
        due:
         dq 2.0
        sei:
         dq 6.0
        fmt:
         db "%lf",10, 0          ; The printf format, "\n",'0'
        msg1: 
        db "Debug",10, 0 

        phi_uno:
        dq -0.00029349
        theta_uno:
        dq 0.995014538


         section .text
        
;Print array routine, first argument in rdi, pointer to array base, lenght in r13
print_array:
        push    rbp
        push    rcx
        
        xor     rcx, rcx                ; rcx fa da contatore, azzeriamolo
        mov     r12,rdi   

print_loop:
        push    rcx                     ; salviamo il registro
        

        ;chiamata a printf
        mov     rdi, fmt                ; formato
        movq    xmm0, [r12+8*rcx]; numero
        mov     rax, 1                  ;abbiamo usato un reistro xmm
        call    printf                  

        pop     rcx                     ; ripristiniamo il registro



        inc     ecx                     ; incrementiamo il contatore
        cmp     rcx,r13                  ; loop sugli elementi
        jb     print_loop                 ; jump se ecx<lenght


        pop     rcx
        pop     rbp
        ret


main:
        push rbx

;inizializzazione del vettore degli n
        xor rcx,rcx        
initialize_n:
        movsd xmm0,[zero]

        cmp rcx,0
        je after
        xor rdx,rdx
mult:
        addsd xmm0, [half]
        inc rdx
        cmp rdx,rcx
        jb mult

after:
        movsd [n_range+8*rcx],xmm0

        inc ecx
        cmp ecx, 8
        jb initialize_n
        
;inizializzazione della mesh
        xor rcx,rcx
initialize_mesh:
        movsd xmm0,[zero]

        cmp rcx,0
        je after_b
        xor rdx,rdx
mult_b:
        addsd xmm0, [h]
        inc rdx
        cmp rdx,rcx
        jb mult_b

after_b: 
        movq [mesh+8*rcx],xmm0

        inc ecx
        cmp ecx, 100
        jb initialize_mesh


;ciclo principale per gli n (che ci vorrà un bel po prima di implementarlo)
main_for:

        ; primi due elementi di theta e phi, determinati a mano per evitare problemi
        xor rax,rax   ;perchè l'ho fatto?
        movq xmm0,[zero]
        movq [phi], xmm0
        movq xmm0,[uno]
        movq [theta], xmm0

        movq xmm0,[phi_uno]
        movq [phi+8], xmm0
        movq xmm0,[theta_uno]
        movq [theta+8], xmm0

        ;lea  rdi,[n_range]
        ;mov r13, 8
        ;call print_array 

        ;lea  rdi,[theta]
        ;mov r13, 2
        ;call print_array

        xor rcx,rcx
        add rcx,2

runge_kutta:
        ;k1-> xmm0 e l1-> xmm1 k2->xmm2 l2->xmm3 k3->xmm4 l3->xmm5 k4->xmm6 l4->xmm7
        mov rax,rcx
        sub rax,1
        ;primo step
        movq xmm0, [phi+8*rax]
        mulsd xmm0,[h]
        movq xmm1, [phi+8*rax]
        mulsd xmm1, [meno_due]
        divsd xmm1, [mesh +8*rax]
        movq xmm2,[theta+8*rax]
        mulsd xmm2,xmm2
        subsd xmm1, xmm2
        mulsd xmm1,[h]
        divsd xmm0,[due]
        divsd xmm1, [due]
        ;secondo step
        xorps xmm2,xmm2
        movq xmm2,[phi+8*rax]
        addsd xmm2,xmm1
        mulsd xmm2, [h]
        movq xmm3,[phi+8*rax]
        addsd xmm3,xmm1
        movq xmm4, [mesh +8*rax] 
        addsd xmm4,[h_half]
        divsd xmm3,xmm4
        mulsd xmm3, [meno_due]
        xorps xmm4,xmm4
        movq xmm4,[theta+8*rax]
        addsd xmm4,xmm0
        mulsd xmm4,xmm4
        subsd xmm3,xmm4
        mulsd xmm3,[h]
        divsd xmm2,[due]
        divsd xmm3, [due]
        ;terzo step
        xorps xmm4,xmm4
        movq xmm4,[phi+8*rax]
        addsd xmm4,xmm3
        mulsd xmm4, [h]
        movq xmm5,[phi+8*rax]
        addsd xmm5,xmm3
        movq xmm6, [mesh +8*rax] 
        addsd xmm6,[h_half]
        divsd xmm5,xmm6
        mulsd xmm5, [meno_due]
        xorps xmm6,xmm6
        movq xmm6,[theta+8*rax]
        addsd xmm6,xmm2
        mulsd xmm6,xmm6
        subsd xmm5,xmm6
        mulsd xmm5,[h]
        ;quarto step
        xorps xmm6,xmm6
        movq xmm6,[phi+8*rax]
        addsd xmm6,xmm5
        mulsd xmm6, [h]
        movq xmm7,[phi+8*rax]
        addsd xmm7,xmm5
        movq xmm8, [mesh +8*rax] 
        addsd xmm8,[h]
        divsd xmm7,xmm8
        mulsd xmm7, [meno_due]
        xorps xmm8,xmm8
        movq xmm8,[theta+8*rax]
        addsd xmm8,xmm4
        mulsd xmm8,xmm8
        subsd xmm7,xmm8
        mulsd xmm7,[h]

        ;abbiamo trovato i nuovi theta e phi! salviamoli nell'array
        xorps xmm8,xmm8
        mulsd xmm2,[due]
        mulsd xmm4,[due]
        movq xmm8,xmm0
        addsd xmm8,xmm2
        addsd xmm8,xmm4
        addsd xmm8,xmm6
        mulsd xmm8,[sei]
        addsd xmm8,[theta+8*rax]
        movq [theta+8*rcx], xmm8

        xorps xmm8,xmm8
        mulsd xmm3,[due]
        mulsd xmm5,[due]
        movq xmm8,xmm1
        addsd xmm8,xmm3
        addsd xmm8,xmm5
        addsd xmm8,xmm7
        mulsd xmm8,[sei]
        addsd xmm8,[phi+8*rax]
        movq [phi+8*rcx], xmm8

        inc rcx
        cmp rcx,100
        jb runge_kutta


        lea  rdi,[theta]
        mov r13, 100
        call print_array

        pop     rbx                     ; ripristino rbx
        ret
       



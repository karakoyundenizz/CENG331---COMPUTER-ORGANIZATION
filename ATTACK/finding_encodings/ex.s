.section __TEXT,__text  // Kodu "text" bölümüne koy
.globl _main            // _main sembolünü global yap (linker için)
_main:                  // Fonksiyonun başlangıç etiketi
    
    // SİZİN KODUNUZ BURADA
    leaq 24(%rsp), %rdi
    leaq 8(%rsp), %rsi

    // Fonksiyondan çıkış
    retq


# modified code from http://rosetta.ahmedmoustafa.io/selection/

diploid_haploid_deriva <- function(selection_coeficient, number_of_generations, number_of_repetitions, freq_recessive, pop_size){
    fitness_aa = 1 - selection_coefficient
    fitness_ab = 1 - selection_coefficient
    fitness_bb = 1
    
    fitness_a = 1 - selection_coefficient
    fitness_b = 1
    
    pdf("simulations_deriva.pdf",)
    par(mar=c(5,5,3,2)+0.1)
    plot(fitness_aa, type="n", xlim=c(0,number_of_generations), ylim=c((freq_recessive-0.1),1), xlab="Generations", ylab="Recessive allele frequency", main="Selection of an allele in Haploid and Diploid populations", cex.sub=1.5, cex.axis=1.5, cex.main=1.5, cex.lab=1.5)
    legend("bottomright", inset=-0.008, legend=c("Haploid Selection", "Diploid Selection"), col=c("Grey 7","Grey 42"), lwd=3, lty=c(4, 1), bty="n", cex=1.5, horiz = T)
        text(c("Selection \n coefficient", selection_coefficient), x=c(number_of_generations/10, number_of_generations/10), y=c(0.98, 0.88), cex=1.5)
    
    for (i in 1 : number_of_repetitions){
        freq_b = freq_recessive
        freq_a = 1 - freq_b
        freq_a_haploid = freq_a
        freq_a_diploid = freq_a
        freq_b_haploid = freq_b
        freq_b_diploid = freq_b
        freq_aa = freq_a^2
        freq_ab = 2 * freq_a * freq_b
        freq_bb = freq_b^2
        
        genotypes = data.frame(t = 0, cbind(freq_aa, freq_ab, freq_bb))
        alleles_diploid = data.frame(t = 0, cbind(freq_a_diploid, freq_b_diploid))
        alleles_haploid = data.frame(t = 0, cbind(freq_a_haploid, freq_b_haploid))
        
        for (t in 1 : number_of_generations) {
            
            x_freq_a = freq_aa + 0.5*freq_ab
            x_freq_b = freq_bb + 0.5*freq_ab
            
            x_freq_aa = freq_aa * fitness_aa
            x_freq_ab = freq_ab * fitness_ab
            x_freq_bb = freq_bb * fitness_bb
            
            sum = x_freq_aa + x_freq_ab + x_freq_bb
            
            freq_aa = x_freq_aa / sum
            freq_ab = x_freq_ab / sum
            freq_bb = x_freq_bb / sum
            
            freq_a_diploid = freq_aa + 0.5 * freq_ab
            freq_b_diploid = freq_bb + 0.5 * freq_ab
            
            A1=rbinom(3,2*pop_size,c(x_freq_bb, x_freq_ab, x_freq_aa))
            x_freq_bb=A1[1]/(pop_size*2)
            x_freq_ab=A1[2]/(pop_size*2)
            x_freq_aa=A1[3]/(pop_size*2)
            
            sum = x_freq_aa + x_freq_ab + x_freq_bb
            
            freq_aa = x_freq_aa / sum
            freq_ab = x_freq_ab / sum
            freq_bb = x_freq_bb / sum
            
            genotypes = rbind(genotypes, data.frame(t, cbind(freq_aa, freq_ab, freq_bb)))
            
            freq_a_diploid = freq_aa + 0.5 * freq_ab
            freq_b_diploid = freq_bb + 0.5 * freq_ab
            
            alleles_diploid = rbind(alleles_diploid, data.frame(t, cbind(freq_a_diploid, freq_b_diploid)))
        }
        lines(alleles_diploid$freq_b, col="Grey 42", lwd=3, cex=1.5)
        
        for (t in 1 : number_of_generations) {
            
            x_freq_a = freq_a_haploid
            x_freq_b = freq_b_haploid
            
            x_freq_a = freq_a_haploid * fitness_a
            x_freq_b = freq_b_haploid * fitness_b
            
            sum = x_freq_a + x_freq_b
            
            freq_a_haploid = x_freq_a / sum
            freq_b_haploid = x_freq_b / sum
            
            x_freq_a=freq_a_haploid
            x_freq_b=freq_b_haploid
            
            A1=rbinom(1,pop_size,x_freq_b)
            x_freq_b=A1/(pop_size);
            x_freq_a=1-x_freq_b            
            
            sum = x_freq_a + x_freq_b
            
            freq_a_haploid = x_freq_a / sum
            freq_b_haploid = x_freq_b / sum
            
            alleles_haploid = rbind(alleles_haploid, data.frame(t, cbind(freq_a_haploid, freq_b_haploid)))
        }
        lines(alleles_haploid$freq_b, col="Grey 7", lwd=3, lty=4, cex=1.5)
        
    }
    dev.off()
}

diploid_haploid_determinista <- function(selection_coeficient, number_of_generations, freq_recessive){
    fitness_aa = 1 - selection_coefficient
    fitness_ab = 1 - selection_coefficient
    fitness_bb = 1
    
    fitness_a = 1 - selection_coefficient
    fitness_b = 1
    
    pdf("simulations_determinista.pdf",)
    par(mar=c(5,5,3,2)+0.1)
    plot(fitness_aa, type="n", xlim=c(0,number_of_generations), ylim=c((freq_recessive-0.1),1), xlab="Generations", ylab="Recessive allele frequency", main="Selection of an allele in Haploid and Diploid populations", cex.sub=1.5, cex.axis=1.5, cex.main=1.5, cex.lab=1.5)
    legend("bottomright", inset=-0.008, legend=c("Haploid Selection", "Diploid Selection"), col=c("Grey 7","Grey 42"), lwd=3, lty=c(4, 1), bty="n", cex=1.5, horiz = T)
    text(c("Selection \n coefficient", selection_coefficient), x=c(number_of_generations/10, number_of_generations/10), y=c(0.98, 0.88), cex=1.5)
    
    for (i in 1 : number_of_repetitions){
        freq_b = freq_recessive
        freq_a = 1 - freq_b
        freq_a_haploid = freq_a
        freq_a_diploid = freq_a
        freq_b_haploid = freq_b
        freq_b_diploid = freq_b
        freq_aa = freq_a^2
        freq_ab = 2 * freq_a * freq_b
        freq_bb = freq_b^2
        
        genotypes = data.frame(t = 0, cbind(freq_aa, freq_ab, freq_bb))
        alleles_diploid = data.frame(t = 0, cbind(freq_a_diploid, freq_b_diploid))
        alleles_haploid = data.frame(t = 0, cbind(freq_a_haploid, freq_b_haploid))
        
        for (t in 1 : number_of_generations) {
            
            x_freq_a = freq_aa + 0.5*freq_ab
            x_freq_b = freq_bb + 0.5*freq_ab
            
            x_freq_aa = freq_aa * fitness_aa
            x_freq_ab = freq_ab * fitness_ab
            x_freq_bb = freq_bb * fitness_bb
            
            sum = x_freq_aa + x_freq_ab + x_freq_bb
            
            freq_aa = x_freq_aa / sum
            freq_ab = x_freq_ab / sum
            freq_bb = x_freq_bb / sum
            
            genotypes = rbind(genotypes, data.frame(t, cbind(freq_aa, freq_ab, freq_bb)))
            
            freq_a_diploid = freq_aa + 0.5 * freq_ab
            freq_b_diploid = freq_bb + 0.5 * freq_ab
            
            alleles_diploid = rbind(alleles_diploid, data.frame(t, cbind(freq_a_diploid, freq_b_diploid)))
        }
        lines(alleles_diploid$freq_b, col="Grey 42", lwd=3, cex=1.5)
        
        for (t in 1 : number_of_generations) {
            
            x_freq_a = freq_a_haploid
            x_freq_b = freq_b_haploid
            
            x_freq_a = freq_a_haploid * fitness_a
            x_freq_b = freq_b_haploid * fitness_b
            
            sum = x_freq_a + x_freq_b
            
            freq_a_haploid = x_freq_a / sum
            freq_b_haploid = x_freq_b / sum
            
            alleles_haploid = rbind(alleles_haploid, data.frame(t, cbind(freq_a_haploid, freq_b_haploid)))
        }
        lines(alleles_haploid$freq_b, col="Grey 7", lwd=3, lty=4, cex=1.5)
        
    }
    dev.off()
}




selection_coefficient = 0.1
number_of_generations = 5*(10^2)
number_of_repetitions = 10^2
freq_recessive = 10^-6
pop_size=10^6

diploid_haploid_deriva(selection_coeficient, number_of_generations, number_of_repetitions, freq_recessive, pop_size)

diploid_haploid_determinista(selection_coeficient, number_of_generations, freq_recessive)

# устанавливаем BiocManager, если его нет
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# устанавливаем нужные пакеты
install.packages('tidyverse')
BiocManager::install(c("DECIPHER", "Biostrings"))

# импортируем пакеты
library(tidyverse)
library(Biostrings)
library(DECIPHER)

# устанавливаем рабочую директорию
setwd('~/R_BLAST/') # заменить на свою

# читаем fasta файл с аминокислотной последовательностью целевого белка у родственного животного
# у нас не было последовательности для ламы, но была для альпаки. ее и использовали
lag3 <- readAAStringSet('lag3_fragment.fasta') # можно не трогать остальные строчки, а только заменить путь к файлу

# BLAST выдал результаты по этой сборке. скачиваем ее и читаем
chromosome <- readDNAStringSet('sequence.fasta') # заменить на путь к файлу 

# так как последовательность нуклеотидная - нужно ее транслировать.
# но у эукариот есть сплайсинг, поэтому чтобы удалить интроны - транслируем в 3-х рамках
# выскакивающие ошибки игнорируем. это из-за некратности длины последовательности 3
frame1 <- translate(chromosome)
frame2 <- translate(subseq(chromosome, 2))
frame3 <- translate(subseq(chromosome, 3))

# делаем выравнивание с нашим референсным белком и сохраняем консенсус
align_frame1 <- AlignSeqs(c(lag3, frame1)) %>% 
  ConsensusSequence() %>% 
  str_remove_all('[+-]') %>% 
  pairwiseAlignment(frame1, .) %>% 
  aligned() %>% 
  str_remove_all('[-+]')

align_frame2 <- AlignSeqs(c(lag3, frame2)) %>% 
  ConsensusSequence() %>% 
  str_remove_all('[+-]') %>% 
  pairwiseAlignment(frame2, .) %>% 
  aligned() %>% 
  str_remove_all('[-+]')

align_frame3 <- AlignSeqs(c(lag3, frame3)) %>% 
  ConsensusSequence() %>% 
  str_remove_all('[+-]') %>% 
  pairwiseAlignment(frame3, .) %>% 
  aligned() %>% 
  str_remove_all('[-+]')

#соединяем последовательности от 3-х рамок в один набор 
consesus_frames <- c(align_frame1, align_frame2, align_frame3)

# делаем выравнивание на референс
sequences <- pairwiseAlignment(consesus_frames, lag3, type='local') %>% 
  aligned() %>% 
  str_remove_all('[-]')

# открываем его в браузере в человекочитаемом виде
alignment <- AlignSeqs(c(lag3, sequences))
alignment %>% 
  BrowseSeqs(colWidth = 150)

# записываем в файл
alignment[2:4] %>% 
  ConsensusSequence() %>% 
  as.character() %>% 
  str_remove_all('[-]') %>% 
  paste('>predicted_protein\n', ., sep = '') %>% 
  write(file='predicted_protein.fasta')



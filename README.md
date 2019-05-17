# motif_finding
The project is designed to search for motifs (equal or similar frequent patterns) and their position weight matrix from artificial gene sequences.

Algorithm Introduction
The frame of greedy algorithm [1] is used in this project. The algorithm starts by forming a pair of most similar candidate motifs from first two sequences. The selected pair is used to find most similar motifs in the thirds sequence, which will join in the pair to build a trine. The join process is repeated in the following sequences until the end of the dataset.
However, we made some modifications:
1.	The first 100 pair of candidate motifs with highest information content will be chosen as starting pair. Evidence suggests that when generating starting pairs, about 10 pairs may reach maximum possible information content. That is, the two components are exactly same by accident. Hence, for cases where the destination motif is not equal, if we only keep the most similar pair, the destination pair will never be chosen as starting pair. 
2.	In each iteration, if the information content of the sites built by any starting pair reaches the maximum possible value, the script will end earlier, as no further candidate set of sites can reach higher information content.

Pseudo Code:
Input: #iteration: loop_num, #stating pairs: startpairs_num, gene sequences: seqs;
Output: predicted sites: pre_sites, predicted PWM: pre_pwm

count = 0, content_info = 0
While count < loop_num:
Disorder seqs
#startpairs_num most similar starting pairs->starting_pairs
For pair in starting_pairs:
temp_sites<-starting_pairs
for seq in seqs[2:]:
temp_sites.append(most similar motif)
ci <-calculate_content_info(temp_sites)
if ci > content_info:
candidate_sites<- temp_sites
if ci== max possible value:
break while
count+=1
pre_pwm <-calculate_pwm(candidate_sites)
return(candidate_sites, pre_pwm)

clc;
clear;


%% Basic information for Vulpes vulpes (Red fox)
mitochondria_gbk = getgenbank('NC_008434');
%load mitochondria
mitochondria = mitochondria_gbk.Sequence;
bases = basecount(mitochondria)
compBases = basecount(seqrcomplement(mitochondria));
figure
ntdensity(mitochondria)

%% The data of all species
data = {'Vulpes vulpes Red Fox'                'NC_008434'   [3 13]  ;
        'Vulpes lagopus Arctic Fox'             'NC_026529'   [3 13]  ;
        'Canis anthus African golden wolf'              'NC_027956'   [3 13]  ;
        'Martes americana American marten'                   'HM106324' [3 13]  ;
        'Nyctereutes procyonoides Raccoon dog'        'MG256392'   [3 13]  ;
        'Prionailurus viverrinus Fishing cat'      'NC_028305'   [3 13]  ;
        'Talpa europaea European mole'   'NC_002391' [3 13]  ;
        'Felis catus Domestic cat'               'NC_001700' [3 13];
         'Leopardus pardalis Ocelot'           'NC_028315' [3 13];
         %'Lynx pardinus Spanish lynx'          'KX911412'   [3 13];
         'Panthera tigris Tiger'                    'EF551003'   [3 13];
         'Canis lupus familiaris dog'            'KF907307'  [3 13];
         'Crocodylus johnsoni river crocodile'  'NC_015238' [3 13];
        };
    
numSpecies = size(data,1);

%%
% Retrieve the sequence information from the NCBI GenBank database for 
% all of the accession numbers. 
%
for ind = 1:numSpecies
     lentivirus(ind) = getgenbank(data{ind,2});
end

% Extract CDS for the CYTB and COX1 coding regions. Then extract the
% nucleotide sequences using the CDS pointers.
for ind = 1:numSpecies
    temp_seq = lentivirus(ind).Sequence; 
    temp_seq = regexprep(temp_seq,'[nry]','a');
    CDSs = lentivirus(ind).CDS(data{ind,3});
    CYTC(ind).Sequence = temp_seq(CDSs(1).indices(1):CDSs(1).indices(2)); 
    CYTB(ind).Sequence = temp_seq(CDSs(2).indices(1):CDSs(2).indices(2)); 
end


%%
% Convert nucleotide sequences to amino acid sequences using |nt2aa|.
for ind = 1:numSpecies
    aacytb(ind).Sequence = nt2aa(CYTB(ind).Sequence);
    aacytc(ind).Sequence = nt2aa(CYTC(ind).Sequence);
end

%% Calculate the distance and linkage, and then generate the tree.

% Tajima-Nei model, which gives a same result, so commented
% cytcnt = seqpdist(CYTC,'method','Tajima-Nei','Alphabet','NT', 'indel','pair');
% cytcnttree = seqlinkage(cytcnt,'UPGMA',data(:,1))
% plot(cytcnttree,'type','angular');
% title('COX1 coding gene')

%CYTB coding gene tree with Kimura model and UPGMA
cytcnt = seqpdist(CYTC,'method','Kimura','Alphabet','NT', 'indel','pair');
cytcnttree = seqlinkage(cytcnt,'UPGMA',data(:,1))
plot(cytcnttree,'type','angular');
title('COX1 coding gene')

%COX1 protein with Jukes-Cantor model and WPGMA
cytcaa = seqpdist(aacytc,'method','Jukes-Cantor','indel','pair');
cytcaatree = seqlinkage(cytcaa,'WPGMA',data(:,1))
plot(cytcaatree,'type','angular');
title('COX1 protein')

%CYTB coding gene tree with Kimura model and UPGMA
cytbnt = seqpdist(CYTB,'method','Kimura','Alphabet','NT', 'indel','pair');
cytbnttree = seqlinkage(cytbnt,'UPGMA',data(:,1))
plot(cytbnttree,'type','angular');
title('CYTB coding gene')

%COX1 protein tree with Jukes-Cantor model and WPGMA
cytbaa = seqpdist(aacytb,'method','Jukes-Cantor','indel','pair');
cytbaatree = seqlinkage(cytbaa,'WPGMA',data(:,1))
plot(cytbaatree,'type','angular');
title('CYTB protein')

%% Genrate weights for consensus tree
weights = [sum(cytbnt) sum(cytbaa) sum(cytcnt) sum(cytcaa)];
weights = weights / sum(weights);
dist = cytbnt .* weights(1) + cytbaa .* weights(2) + cytcnt .* weights(3) + cytcaa .* weights(4);

%%Generate consensus tree
con_tree = seqlinkage(dist,'average',data(:,1));
plot(con_tree,'type','angular');
title('Consensus weighted tree')
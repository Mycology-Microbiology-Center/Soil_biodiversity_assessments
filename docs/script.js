// Initialize Vue with Vuetify
new Vue({
    el: '#app',
    vuetify: new Vuetify(),
    data: {
        title: 'Experiment Cost Calculator',
        
        // Experimental design data
        poolingMethod: ['Unpooled', 'DNA Pooling', 'Soil pooling', 'Semi-pooling'],
        numSites: 1,
        numSamples: 10,
        numSemipools: 9,
        selectedObject: null,

        
        // Cost inputs
        dnaExtractionCost: 5,
        pcrCost: 3,
        librarySequencingCost: 2000,
        
        // Sequencing parameters
        requiredReads: 10000,
        qualityControlLossRate: 20,
        poolingFactor: 0.5,
        sequencingPlatform: 'PacBio Sequel II',
        platforms: ['Illumina MiSeq', 'Illumina NextSeq', 'Illumina NovaSeq', 'PacBio Sequel II', 'PacBio Revio', 'Oxford Nanopore', 'Element Biosciences AVITI', 'MGI DNBSEQ-G99', 'Other'],
        sequencingThroughput: 8000000,
        platformThroughputs: {
            'Illumina MiSeq':            25000000,    // ~25M reads per run (V3 600 cycle kit)
            'Illumina NextSeq':          400000000,   // ~400M reads per run (high output)
            'Illumina NovaSeq':          10000000000, // ~10B reads per run (S4 flow cell)
            'PacBio Sequel II':          8000000,     // ~8M HiFi reads per SMRT cell
            'PacBio Revio':              25000000,    // ~25M HiFi reads per SMRT cell
            'Oxford Nanopore':           50000000,    // ~50M reads for PromethION flow cell
            'Element Biosciences AVITI': 100000000,   // ~100M reads per run ('Medium Output' flow cell)
            'MGI DNBSEQ-G99':            80000000,    // ~80M reads per run (FCL flow cell type)
            'Other':                     15000000     // Default value for a user-defined platform
        },
        userOverrideThroughput: false,      // Track if user manually changed the throughput value
        
        // Results data
        resultHeaders: [
            { text: 'Metric', value: 'metric', align: 'start', sortable: false },
            { text: 'Unpooled', value: 'unpooled', align: 'center', sortable: false },
            { text: 'DNA Pooling', value: 'dnaPooling', align: 'center', sortable: false },
            { text: 'Soil pooling', value: 'soilPooling', align: 'center', sortable: false },
            { text: 'Semi-pooling', value: 'semiPooling', align: 'center', sortable: false }
                ]
    },

    computed: {
        recommendedPoolingFactor() {
          if (!this.selectedObject) return 1;
          switch (this.selectedObject) {
            case 'Bacteria': return 0.25;
            case 'Fungi': return 0.75;
            case 'Animal': return 1;
            default: return 1;
          }
        },
        
        qualityLossRateDecimal() {
          // Convert percentage to decimal (e.g., 50% -> 0.5)
          return this.qualityControlLossRate / 100;
        },
        
        sequencingCapacity() {
          if (!this.sequencingThroughput || !this.requiredReads) {
            return {
              unpooled: { maxSamples: 0, currentSamples: 0, remainingCapacity: 0, utilizationPercent: 0 },
              pooled: { maxSamples: 0, currentSamples: 0, remainingCapacity: 0, utilizationPercent: 0 },
              semiPooled: { maxSamples: 0, currentSamples: 0, remainingCapacity: 0, utilizationPercent: 0 }
            };
          }
          
          const adjustedReads = this.requiredReads + (this.requiredReads * this.qualityLossRateDecimal);
          
          // Calculate capacity for each method
          const methods = {
            unpooled: {
              sequencingUnits: this.numSites * this.numSamples,
              readsPerUnit: adjustedReads
            },
            pooled: {
              sequencingUnits: this.numSites,
              readsPerUnit: this.numSamples * this.poolingFactor * adjustedReads
            },
            semiPooled: {
              sequencingUnits: this.numSites * this.numSemipools,
              readsPerUnit: this.numSemipools > 0 ? (this.numSamples / this.numSemipools) * this.poolingFactor * adjustedReads : 0
            }
          };
          
          const result = {};
          for (const [method, data] of Object.entries(methods)) {
            const maxSamples = data.readsPerUnit > 0 ? Math.floor(this.sequencingThroughput / data.readsPerUnit) : 0;
            const currentSamples = data.sequencingUnits;
            const remainingCapacity = Math.max(0, maxSamples - currentSamples);
            const utilizationPercent = maxSamples > 0 ? Math.round((currentSamples / maxSamples) * 100) : 0;
            
            result[method] = {
              maxSamples,
              currentSamples,
              remainingCapacity,
              utilizationPercent
            };
          }
          
          return result;
        },
        
        isValidConfiguration() {
          return this.sequencingPlatform && 
                 this.numSites >= 1 && 
                 this.numSamples >= 1 &&
                 this.numSemipools >= 1 &&
                 this.requiredReads > 0 &&
                 this.qualityControlLossRate >= 0 &&
                 this.qualityControlLossRate < 100 &&
                 this.sequencingThroughput > 0 &&
                 this.dnaExtractionCost >= 0 &&
                 this.pcrCost >= 0 &&
                 this.librarySequencingCost >= 0;
        },
        
        totalSamples() {
          return this.numSites * this.numSamples;
        },
        
        // Distribute samples into semipools
        semipoolBalance() {
          if (this.numSemipools <= 0 || this.numSamples <= 0) {
            return { isBalanced: true, message: '' };
          }
          
          const samplesPerSemipool = Math.floor(this.numSamples / this.numSemipools);
          const remainder = this.numSamples % this.numSemipools;
          
          if (remainder === 0) {
            return { 
              isBalanced: true, 
              message: `${samplesPerSemipool} samples per semipool` 
            };
          } else {
            return { 
              isBalanced: false, 
              message: `Uneven distribution: ${remainder} semipool(s) will have ${samplesPerSemipool + 1} samples, ${this.numSemipools - remainder} will have ${samplesPerSemipool} samples` 
            };
          }
        },
        
        calculatedResults() {
          if (!this.isValidConfiguration) return [];
          
          const results = {
            totalReads: {
              metric: 'Total number of reads required',
              unpooled: this.calculateTotalReads('Unpooled', this.totalSamples),
              dnaPooling: this.calculateTotalReads('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateTotalReads('Soil pooling', this.totalSamples),
              semiPooling: this.calculateTotalReads('Semi-pooling', this.totalSamples),
            },
            sequencingRuns: {
              metric: 'Number of sequencing runs',
              unpooled: this.calculateSequencingRuns('Unpooled', this.totalSamples),
              dnaPooling: this.calculateSequencingRuns('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateSequencingRuns('Soil pooling', this.totalSamples),
              semiPooling: this.calculateSequencingRuns('Semi-pooling', this.totalSamples),
            },
            total: {
              metric: 'Total cost:',
              unpooled: this.calculateTotalCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculateTotalCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateTotalCost('Soil pooling', this.totalSamples),
              semiPooling: this.calculateTotalCost('Semi-pooling', this.totalSamples),
            },
            costExtraction: {
              metric: '    • DNA extraction',
              unpooled: this.calculateDnaExtractionCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculateDnaExtractionCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateDnaExtractionCost('Soil pooling', this.totalSamples),
              semiPooling: this.calculateDnaExtractionCost('Semi-pooling', this.totalSamples),
            },
            costPCR: {
              metric: '    • PCR',
              unpooled: this.calculatePcrCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculatePcrCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculatePcrCost('Soil pooling', this.totalSamples),
              semiPooling: this.calculatePcrCost('Semi-pooling', this.totalSamples),
            },
            costSequencing: {
              metric: '    • Sequencing',
              unpooled: this.calculateSequencingCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculateSequencingCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateSequencingCost('Soil pooling', this.totalSamples),
              semiPooling: this.calculateSequencingCost('Semi-pooling', this.totalSamples),
            }
          };
          
          return [
            results.totalReads,
            results.sequencingRuns,
            results.total,
            results.costExtraction,
            results.costPCR,
            results.costSequencing
          ];
        }
      },
      
    
    methods: {
        onObjectChange() {     
            if (this.selectedObject === 'Bacteria') {
            this.poolingFactor = 0.25;
            } else if (this.selectedObject === 'Fungi') {
            this.poolingFactor = 0.75;
            } else if (this.selectedObject === 'Animal') {
            this.poolingFactor = 1;
            }
                },

        // Total number of reads required
        calculateTotalReads(method, totalSamples) {
            // Account for quality control losses by adjusting required reads
            const adjustedReads = this.requiredReads + (this.requiredReads * this.qualityLossRateDecimal);
            
            switch(method) {
                case 'Unpooled':
                    // Each sample sequenced individually
                    return Math.ceil(totalSamples * adjustedReads);
                    
                case 'DNA Pooling':
                    // DNA from all samples within each site is pooled into one pooled sample per site
                    // Total reads = sites × (samples per site × pooling factor × adjusted reads per sample)
                    return Math.ceil(this.numSites * (this.numSamples * this.poolingFactor * adjustedReads));
                    
                case 'Soil pooling':
                    // Soil samples are pooled before DNA extraction within each site
                    // Total reads = sites × (samples per site × pooling factor × adjusted reads per sample)
                    return Math.ceil(this.numSites * (this.numSamples * this.poolingFactor * adjustedReads));
                    
                case 'Semi-pooling':
                    // Samples within each site are partitioned into semipools, each semipool processed separately
                    // Total reads = sites × semipools per site × (samples per semipool × pooling factor × adjusted reads per sample)
                    // This simplifies to: sites × numSamples × pooling factor × adjusted reads per sample
                    return Math.ceil(this.numSites * (this.numSamples * this.poolingFactor * adjustedReads));
                    
                default:
                    return 0;
            }
        },

        // Number of sequencing runs required
        calculateSequencingRuns(method, totalSamples) {
            const totalReads = this.calculateTotalReads(method, totalSamples);
            const runs = Math.ceil(totalReads / this.sequencingThroughput);
            return runs;
        },
        
        // DNA extraction cost
        calculateDnaExtractionCost(method, totalSamples) {
            let dnaExtractionSamples;
            
            switch(method) {
                case 'Unpooled':
                    // Extract DNA from each individual sample
                    dnaExtractionSamples = totalSamples;
                    break;
                    
                case 'DNA Pooling':
                    // Extract DNA from each individual sample (pooling happens after extraction)
                    dnaExtractionSamples = totalSamples;
                    break;
                    
                case 'Soil pooling':
                    // Pool soil samples before extraction, so only extract from pooled samples (one per site)
                    dnaExtractionSamples = this.numSites;
                    break;
                    
                case 'Semi-pooling':
                    // Samples are pooled into semipools before extraction (one extraction per semipool)
                    dnaExtractionSamples = this.numSites * this.numSemipools;
                    break;
                    
                default:
                    dnaExtractionSamples = 0;
            }
            
            return dnaExtractionSamples * this.dnaExtractionCost;
        },
          
        // PCR cost
        calculatePcrCost(method, totalSamples) {
            let pcrSamples;
            
            switch(method) {
                case 'Unpooled':
                    // PCR each individual sample
                    pcrSamples = totalSamples;
                    break;
                    
                case 'DNA Pooling':
                    // Pool DNA within each site, then PCR the pooled samples (one PCR per site)
                    pcrSamples = this.numSites;
                    break;
                    
                case 'Soil pooling':
                    // Soil is pooled before extraction, so one PCR per site
                    pcrSamples = this.numSites;
                    break;
                    
                case 'Semi-pooling':
                    // Each semipool is PCR'd separately (one PCR per semipool)
                    pcrSamples = this.numSites * this.numSemipools;
                    break;
                    
                default:
                    pcrSamples = 0;
            }
            
            return pcrSamples * this.pcrCost;
        },
          
        // Sequencing cost
        calculateSequencingCost(method, totalSamples) {
            // Sequencing cost is based on the number of sequencing runs needed
            const sequencingRuns = this.calculateSequencingRuns(method, totalSamples);
            return sequencingRuns * this.librarySequencingCost;
        },
          
        // Total cost
        calculateTotalCost(method, totalSamples) {
            // Calculate total cost by summing all component costs
            const dnaExtractionCost = this.calculateDnaExtractionCost(method, totalSamples);
            const pcrCost = this.calculatePcrCost(method, totalSamples);
            const sequencingCost = this.calculateSequencingCost(method, totalSamples);
            
            return Math.round(dnaExtractionCost + pcrCost + sequencingCost);
        },
        
        updateThroughput() {
            // Only update if user hasn't manually changed the value
            if (!this.userOverrideThroughput && this.sequencingPlatform) {
                this.sequencingThroughput = this.platformThroughputs[this.sequencingPlatform];
            }
        },
        
        onThroughputInput() {
            // Mark that user has manually changed the throughput
            this.userOverrideThroughput = true;
        },
        
        resetThroughputOverride() {
            // Reset the override flag and update with platform default
            this.userOverrideThroughput = false;
            this.updateThroughput();
        }
    },
    watch: {
        sequencingPlatform(newValue) {
            if (newValue) {
                this.updateThroughput();
            }
        }
    },
    mounted() {
        console.log('app mounted successfully!');
    }
}); 

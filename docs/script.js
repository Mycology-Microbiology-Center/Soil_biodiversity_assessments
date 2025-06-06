// Initialize Vue with Vuetify
new Vue({
    el: '#app',
    vuetify: new Vuetify(),
    data: {
        title: 'Experiment Cost Calculator',
        
        // Experimental design data
        poolingMethod: ['Unpooled', 'DNA Pooling', 'Soil pooling'],
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
        poolingEffect: 0.5,
        sequencingPlatform: 'PacBio Sequel II',
        platforms: ['Illumina MiSeq', 'Illumina NextSeq', 'Illumina NovaSeq', 'PacBio Sequel II', 'PacBio Revio', 'Oxford Nanopore', 'Other'],
        sequencingThroughput: 8000000,
        platformThroughputs: {
            'Illumina MiSeq':   25000000,    // ~25M reads per run (V3 600 cycle kit)
            'Illumina NextSeq': 400000000,   // ~400M reads per run (high output)
            'Illumina NovaSeq': 10000000000, // ~10B reads per run (S4 flow cell)
            'PacBio Sequel II': 8000000,     // ~8M HiFi reads per SMRT cell
            'PacBio Revio':     25000000,    // ~25M HiFi reads per SMRT cell
            'Oxford Nanopore':  50000000,    // ~50M reads for PromethION flow cell
            'Other':            15000000     // Default value for a user-defined platform
        },
        userOverrideThroughput: false,      // Track if user manually changed the throughput value
        
        // Results data
        resultHeaders: [
            { text: 'Metric', value: 'metric', align: 'start', sortable: false },
            { text: 'Unpooled', value: 'unpooled', align: 'center', sortable: false },
            { text: 'DNA Pooling', value: 'dnaPooling', align: 'center', sortable: false },
            { text: 'Soil pooling', value: 'soilPooling', align: 'center', sortable: false }
                ]
    },

    computed: {
        recommendedPoolingEffect() {
          if (!this.selectedObject) return 1;
          switch (this.selectedObject) {
            case 'Bacteria': return 0.25;
            case 'Fungi': return 0.75;
            case 'Animal': return 1;
            default: return 1;
          }
        },
        
        sequencingCapacity() {
          if (!this.sequencingThroughput || !this.requiredReads) {
            return {
              maxSamples: 0,
              currentSamples: 0,
              remainingCapacity: 0,
              utilizationPercent: 0
            };
          }
          
          // Calculate initial reads needed accounting for quality filtering losses
          // Assuming ~50% of reads may be lost to quality filtering, chimeras, non-specific amplification
          const qualityLossRate = 0.5; // 50% loss rate
          const initialReadsNeeded = this.requiredReads / (1 - qualityLossRate);
          const maxSamples = Math.floor(this.sequencingThroughput / initialReadsNeeded);
          const currentSamples = this.numSites * this.numSamples;
          const remainingCapacity = Math.max(0, maxSamples - currentSamples);
          const utilizationPercent = currentSamples > 0 ? Math.round((currentSamples / maxSamples) * 100) : 0;
          
          return {
            maxSamples,
            currentSamples,
            remainingCapacity,
            utilizationPercent
          };
        },
        
        isValidConfiguration() {
          return this.sequencingPlatform && 
                 this.numSites >= 1 && 
                 this.numSamples >= 1 &&
                 this.numSemipools >= 1 &&
                 this.requiredReads > 0 &&
                 this.sequencingThroughput > 0 &&
                 this.dnaExtractionCost >= 0 &&
                 this.pcrCost >= 0 &&
                 this.librarySequencingCost >= 0;
        },
        
        totalSamples() {
          return this.numSites * this.numSamples;
        },
        
        calculatedResults() {
          if (!this.isValidConfiguration) return [];
          
          const results = {
            totalReads: {
              metric: 'Total number of reads required',
              unpooled: this.calculateTotalReads('Unpooled', this.totalSamples),
              dnaPooling: this.calculateTotalReads('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateTotalReads('Soil pooling', this.totalSamples),
            },
            sequencingRuns: {
              metric: 'Number of sequencing runs',
              unpooled: this.calculateSequencingRuns('Unpooled', this.totalSamples),
              dnaPooling: this.calculateSequencingRuns('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateSequencingRuns('Soil pooling', this.totalSamples),
            },
            total: {
              metric: 'Total cost:',
              unpooled: this.calculateTotalCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculateTotalCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateTotalCost('Soil pooling', this.totalSamples),
            },
            costExtraction: {
              metric: '    • DNA extraction',
              unpooled: this.calculateDnaExtractionCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculateDnaExtractionCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateDnaExtractionCost('Soil pooling', this.totalSamples),
            },
            costPCR: {
              metric: '    • PCR',
              unpooled: this.calculatePcrCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculatePcrCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculatePcrCost('Soil pooling', this.totalSamples),
            },
            costSequencing: {
              metric: '    • Sequencing',
              unpooled: this.calculateSequencingCost('Unpooled', this.totalSamples),
              dnaPooling: this.calculateSequencingCost('DNA Pooling', this.totalSamples),
              soilPooling: this.calculateSequencingCost('Soil pooling', this.totalSamples),
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
            this.poolingEffect = 0.25;
            } else if (this.selectedObject === 'Fungi') {
            this.poolingEffect = 0.75;
            } else if (this.selectedObject === 'Animal') {
            this.poolingEffect = 1;
            }
                },

        calculateTotalReads(method, totalSamples) {
            let samples, poolFactor;
            
            switch(method) {
                case 'Unpooled':
                    return Math.ceil(totalSamples * this.requiredReads);
                case 'DNA Pooling':
                    samples = totalSamples;
                    return Math.ceil(this.requiredReads * (this.poolingEffect * this.numSemipools));
                case 'Soil pooling':
                    samples = this.numSites;
                    return Math.ceil(this.requiredReads * (this.poolingEffect * this.numSemipools));
                 default:
                    return 0;
            }

        },



        calculateSequencingRuns(method, totalSamples) {
            const totalReads = this.calculateTotalReads(method, totalSamples);
            const runs = Math.ceil(totalReads / this.sequencingThroughput);
            return runs;
        },
        
        calculateDnaExtractionCost(method, totalSamples) {
            let dnaExtractionSamples, pcrSamples, sequencingRuns;
            
            switch(method) {
                case 'Unpooled':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = totalSamples;
                    break;
                case 'DNA Pooling':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = Math.max(1, Math.floor(totalSamples / this.numSemipools));
                    break;
                case 'Soil pooling':
                    dnaExtractionSamples = this.numSites;
                    pcrSamples = this.numSites;
                    break;
                default:
                    dnaExtractionSamples = 0;
                    pcrSamples = 0;
            }
            return dnaExtractionSamples * this.dnaExtractionCost;
          },
          
          calculatePcrCost(method, totalSamples) {
            let dnaExtractionSamples, pcrSamples, sequencingRuns;
            
            switch(method) {
                case 'Unpooled':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = totalSamples;
                    break;
                case 'DNA Pooling':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = Math.max(1, Math.floor(totalSamples / this.numSemipools));
                    break;
                case 'Soil pooling':
                    dnaExtractionSamples = this.numSites;
                    pcrSamples = this.numSites;
                    break;
                default:
                    dnaExtractionSamples = 0;
                    pcrSamples = 0;
            }
            return pcrSamples * this.pcrCost;
          },
          
          calculateSequencingCost(method, totalSamples) {
            let dnaExtractionSamples, pcrSamples, sequencingRuns;
            
            switch(method) {
                case 'Unpooled':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = totalSamples;
                    break;
                case 'DNA Pooling':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = Math.max(1, Math.floor(totalSamples / this.numSemipools));
                    break;
                case 'Soil pooling':
                    dnaExtractionSamples = this.numSites;
                    pcrSamples = this.numSites;
                    break;
                default:
                    dnaExtractionSamples = 0;
                    pcrSamples = 0;
            }
            sequencingRuns = this.calculateSequencingRuns(method, totalSamples);
            return sequencingRuns * this.librarySequencingCost;
          },
          
        calculateTotalCost(method, totalSamples) {
            let dnaExtractionSamples, pcrSamples, sequencingRuns;
            
            switch(method) {
                case 'Unpooled':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = totalSamples;
                    break;
                case 'DNA Pooling':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = Math.max(1, Math.floor(totalSamples / this.numSemipools));
                    break;
                case 'Soil pooling':
                    dnaExtractionSamples = this.numSites;
                    pcrSamples = this.numSites;
                    break;
                default:
                    dnaExtractionSamples = 0;
                    pcrSamples = 0;
            }
            
            sequencingRuns = this.calculateSequencingRuns(method, totalSamples);
            
            const dnaExtractCost = dnaExtractionSamples * this.dnaExtractionCost;
            const pcrCost = pcrSamples * this.pcrCost;
            const sequencingCost = sequencingRuns * this.librarySequencingCost;
            
            return  Math.round(dnaExtractCost + pcrCost + sequencingCost);
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

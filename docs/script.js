// Initialize Vue with Vuetify
new Vue({
    el: '#app',
    vuetify: new Vuetify(),
    data: {
        title: 'Experiment Cost Calculator',
        message: 'Configure your parameters and click Calculate',
        isBlue: true,
        
        // Experimental design data
        poolingMethod: ['Unpooled', 'DNA Pooling', 'Soil pooling', 'Semi-pooled'],
        numSites: 1,
        numSamples: 10,
        numSemipools: 9,
        selectedObject: 'Fungi',

        
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
        showResults: false,
        resultHeaders: [
            { text: 'Metric', value: 'metric', align: 'start', sortable: false },
            { text: 'Unpooled', value: 'unpooled', align: 'center', sortable: false },
            { text: 'DNA Pooling', value: 'dnaPooling', align: 'center', sortable: false },
            { text: 'Soil pooling', value: 'soilPooling', align: 'center', sortable: false },
            { text: 'Semi-pooled', value: 'semiPooled', align: 'center', sortable: false }
        ],
        resultItems: []
    },
    methods: {
        buttonClicked() {
            // Validate form
            if (!this.sequencingPlatform) {
                this.message = 'Please select a sequencing platform';
                this.isBlue = false;
                return;
            }
            
            const totalSamples = this.numSites * this.numSamples;
            
            // Calculate results for each pooling method
            const results = {
                totalReads: {
                    metric: 'Total number of reads required',
                    unpooled: this.calculateTotalReads('Unpooled', totalSamples),
                    dnaPooling: this.calculateTotalReads('DNA Pooling', totalSamples),
                    soilPooling: this.calculateTotalReads('Soil pooling', totalSamples),
                    semiPooled: this.calculateTotalReads('Semi-pooled', totalSamples)
                },
                sequencingRuns: {
                    metric: 'Number of sequencing runs',
                    unpooled: this.calculateSequencingRuns('Unpooled', totalSamples),
                    dnaPooling: this.calculateSequencingRuns('DNA Pooling', totalSamples),
                    soilPooling: this.calculateSequencingRuns('Soil pooling', totalSamples),
                    semiPooled: this.calculateSequencingRuns('Semi-pooled', totalSamples)
                },
                cost: {
                    metric: 'Expected cost (€)',
                    unpooled: this.calculateTotalCost('Unpooled', totalSamples),
                    dnaPooling: this.calculateTotalCost('DNA Pooling', totalSamples),
                    soilPooling: this.calculateTotalCost('Soil pooling', totalSamples),
                    semiPooled: this.calculateTotalCost('Semi-pooled', totalSamples)
                }
            };
            
            this.resultItems = [
                results.totalReads,
                results.sequencingRuns,
                results.cost
            ];
            
            this.showResults = true;
            this.message = 'Calculation completed successfully. See comparison table below.';
            this.isBlue = true;
        },
        
        calculateTotalReads(method, totalSamples) {
            let samples, poolFactor;
            
            switch(method) {
                case 'Unpooled':
                    return totalSamples * this.requiredReads;
                case 'DNA Pooling':
                    samples = totalSamples;
                    poolFactor = Math.max(1, Math.floor(samples / this.numSemipools));
                    return this.requiredReads * (this.poolingEffect * samples);
                case 'Soil pooling':
                    samples = this.numSites;
                    return this.requiredReads * (this.poolingEffect * totalSamples);
                case 'Semi-pooled':
                    samples = this.numSites * this.numSemipools;
                    return this.requiredReads * (this.poolingEffect * totalSamples);
                default:
                    return 0;
            }

        },

        onObjectChange() {     
                    if (this.selectedObject === 'Bacteria') {
                    this.poolingEffect = 0.25;
                    } else if (this.selectedObject === 'Fungi') {
                    this.poolingEffect = 0.75;
                    } else if (this.selectedObject === 'Animal') {
                    this.poolingEffect = 1;
                    }
                        }
        
        calculateSequencingRuns(method, totalSamples) {
            const totalReads = this.calculateTotalReads(method, totalSamples);
            const runs = Math.ceil(totalReads / this.sequencingThroughput);
            return runs;
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
                case 'Semi-pooled':
                    dnaExtractionSamples = totalSamples;
                    pcrSamples = this.numSites * this.numSemipools;
                    break;
                default:
                    dnaExtractionSamples = 0;
                    pcrSamples = 0;
            }
            
            sequencingRuns = this.calculateSequencingRuns(method, totalSamples);
            
            const dnaExtractCost = dnaExtractionSamples * this.dnaExtractionCost;
            const pcrCost = pcrSamples * this.pcrCost;
            const sequencingCost = sequencingRuns * this.librarySequencingCost;
            
            return Math.round(dnaExtractCost + pcrCost + sequencingCost);
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

// Initialize Vue with Vuetify
new Vue({
    el: '#app',
    vuetify: new Vuetify(),
    data: {
        title: 'Experiment Cost Calculator',
        message: 'Configure your parameters and click Calculate',
        isBlue: true,
        
        // Form data
        selectedStrategy: null,
        poolingMethod: ['Soil pooling', 'Unpooled', 'DNA Pooling'],
        organism: null,
        organisms: ['Bacteria', 'Fungi', 'Animal'],
        numSites: 1,
        numSamples: 10,
        pcrTime: 1,
        samplesPerPCRrun: 96,
        dnaExtractionTime: 2,
        samplesPerDNAextraction: 24,
        personalHourRate: 0,
        poolingNumber: 9,
        samplesPerLibrary: 1,
        costOfLibrary: null,
        pricePerSample: null,
        sequencingDepth: 100,
        sampleSequencingDepth: 10,
        x: 10,
        // Consumables table
        consumableHeaders: [
            { text: 'Item', value: 'name', sortable: true },
            { text: 'Category', value: 'category', sortable: true },
            { text: 'Price ($)', value: 'price', sortable: true },
            { text: 'Total Volume', value: 'totalVolume', sortable: true },
            { text: 'Volume per Sample', value: 'volumePerSample', sortable: true },
            { text: '', value: 'actions', sortable: false }
        ],
        consumables: [],
        consumableDialog: false,
        editedIndex: -1,
        editedItem: {
            name: '',
            category: '',
            price: 0,
            totalVolume: 0,
            volumePerSample: 0
        },
        defaultItem: {
            name: '',
            category: '',
            price: 0,
            totalVolume: 0,
            volumePerSample: 0
        },
        consumableWarnings: [],
        showConsumableWarningDialog: false
    },
    computed: {
        formTitle() {
            return this.editedIndex === -1 ? 'New Consumable' : 'Edit Consumable'
        },
        sequencingDepthConfig() {
            if (!this.selectedStrategy || !this.organism) return { value: 100, range: false, min: 25, max: 100, recommendedMax: 100, disabled: true };
            
            if (this.selectedStrategy === 'Unpooled') {
                return { value: 100, range: false, min: 25, max: 100, recommendedMax: 100, disabled: true };
            }
            
            switch(this.organism) {
                case 'Animal':
                    return { value: 100, range: false, min: 25, max: 100, recommendedMax: 100 };
                case 'Bacteria':
                    return { value: 25, range: true, min: 25, max: 100, recommendedMax: 75, recommendedMin: 25 };
                case 'Fungi':
                    if (this.selectedStrategy === 'Soil pooling') {
                        return { value: 100, range: false, min: 25, max: 100, recommendedMax: 100 };
                    } else if (this.selectedStrategy === 'DNA Pooling') {
                        return { value: 50, range: true, min: 25, max: 100, recommendedMax: 100, recommendedMin: 50 };
                    }
                    return { value: 100, range: false, min: 25, max: 100, recommendedMax: 100 };
                default:
                    return { value: 100, range: false, min: 25, max: 100, recommendedMax: 100 };
            }
        }
    },
    methods: {
        buttonClicked() {
            // Validate form
            if (!this.selectedStrategy) {
                this.message = 'Please select a strategy';
                this.isBlue = false;
                return;
            }
            
           let nrOfSamplesForExtraction
           let nrOfSamplesForPCR

           if (this.selectedStrategy === 'Soil pooling') {
            nrOfSamplesForExtraction = this.numSites * this.numSamples / this.poolingNumber;
           } else {
            nrOfSamplesForExtraction = this.numSites * this.numSamples;
           }

           if (this.selectedStrategy === 'Unpooled') {
            nrOfSamplesForPCR = this.numSites * this.numSamples;
           } else {
            nrOfSamplesForPCR = this.numSites * this.numSamples / this.poolingNumber;
           }    
            
           let dnaExtractionHours = Math.ceil(nrOfSamplesForExtraction / this.samplesPerDNAextraction) * this.dnaExtractionTime;
           let pcrHours = Math.ceil(nrOfSamplesForPCR / this.samplesPerPCRrun) * this.pcrTime;
           let laborCost = (pcrHours + dnaExtractionHours) * this.personalHourRate;
           let consumablesCost = this.calculateConsumablesCost(nrOfSamplesForPCR, nrOfSamplesForExtraction);
           let libraryCost = this.costOfLibrary * Math.ceil((nrOfSamplesForPCR / this.samplesPerLibrary));
           let sequencingCost = this.calculateSequencingCost(nrOfSamplesForPCR);    
           let totalCost = laborCost + consumablesCost + sequencingCost;
            
            // Show warning dialog if any consumable is insufficient
            if (this.consumableWarnings.length > 0) {
                this.showConsumableWarningDialog = true;
            }
            
            this.message = `Total cost estimate: $${totalCost}\n` +
                          `Labor cost: $${laborCost} (${pcrHours + dnaExtractionHours} hours)\n` +
                          `- DNA extraction: $${dnaExtractionHours * this.personalHourRate} (${dnaExtractionHours} hours)\n` +
                          `- PCR work: $${pcrHours * this.personalHourRate} (${pcrHours} hours)\n` +
                          `Consumables cost: $${consumablesCost} (${nrOfSamplesForPCR} samples)\n` +
                          `Library cost: $${libraryCost.toFixed(1)}\n` +
                          `Sequencing cost: $${sequencingCost.toFixed(1)} (${nrOfSamplesForPCR} samples)`;
            this.isBlue = true;
            
            console.log('Form submitted:', {
                strategy: this.selectedStrategy,
                sites: this.numSites,
                samples: this.numSamples,
                nrOfSamplesForExtraction: nrOfSamplesForExtraction,
                nrOfSamplesForPCR: nrOfSamplesForPCR,
                dnaExtractionHours: dnaExtractionHours,
                pcrHours: pcrHours,
                laborCost: laborCost,
                consumablesCost: this.calculateConsumablesCost(nrOfSamplesForPCR),
                totalCost: totalCost
            });
        },

        // Consumables table methods
        addConsumable() {
            this.editedIndex = -1;
            this.editedItem = Object.assign({}, this.defaultItem);
            this.consumableDialog = true;
        },

        editConsumable(item) {
            this.editedIndex = this.consumables.indexOf(item);
            this.editedItem = Object.assign({}, item);
            this.consumableDialog = true;
        },

        deleteConsumable(item) {
            const index = this.consumables.indexOf(item);
            if (confirm('Are you sure you want to delete this item?')) {
                this.consumables.splice(index, 1);
            }
        },

        closeConsumableDialog() {
            this.consumableDialog = false;
            this.$nextTick(() => {
                this.editedItem = Object.assign({}, this.defaultItem);
                this.editedIndex = -1;
            });
        },

        saveConsumable() {
            if (!this.editedItem.name || !this.editedItem.category || !this.editedItem.price || !this.editedItem.totalVolume || !this.editedItem.volumePerSample) {
                alert('Please fill in all fields, including category');
                return;
            }

            if (this.editedIndex > -1) {
                Object.assign(this.consumables[this.editedIndex], this.editedItem);
            } else {
                this.consumables.push(this.editedItem);
            }
            this.closeConsumableDialog();
        },

        calculateConsumablesCost(nrOfSamplesForPCR, nrOfSamplesForExtraction) {
            this.consumableWarnings = [];
            return this.consumables.reduce((total, item) => {
                let relevantSamples = 0;
                if (item.category === 'DNA extraction') {
                    relevantSamples = nrOfSamplesForExtraction;
                } else if (item.category === 'PCR') {
                    relevantSamples = nrOfSamplesForPCR;
                }
                const totalVolumeNeeded = item.volumePerSample * relevantSamples;
                if (totalVolumeNeeded > item.totalVolume) {
                    this.consumableWarnings.push(`Not enough '${item.name}' for ${item.category}: need ${totalVolumeNeeded}, have ${item.totalVolume}`);
                }
                const unitsNeeded = totalVolumeNeeded / item.totalVolume;
                return total + (unitsNeeded * item.price);
            }, 0);
        },
        calculateSequencingCost(nrOfSamplesForPCR) {
            let sequencingCost = nrOfSamplesForPCR * this.pricePerSample;
            return sequencingCost;
        }
    },
    watch: {
        organism(newValue) {
            if (this.selectedStrategy === 'Unpooled') {
                this.sequencingDepth = 100;
                return;
            }
            
            switch(newValue) {
                case 'Bacteria':
                    this.sequencingDepth = 25;
                    break;
                case 'Animal':
                    this.sequencingDepth = 100;
                    break;
                case 'Fungi':
                    if (this.selectedStrategy === 'Soil pooling') {
                        this.sequencingDepth = 100;
                    } else if (this.selectedStrategy === 'DNA Pooling') {
                        this.sequencingDepth = 50;
                    } else {
                        this.sequencingDepth = 100;
                    }
                    break;
                default:
                    this.sequencingDepth = 100;
            }
        },
        selectedStrategy(newValue) {
            if (newValue === 'Unpooled') {
                this.sequencingDepth = 100;
                return;
            }
            
            if (this.organism === 'Fungi') {
                if (newValue === 'Soil pooling') {
                    this.sequencingDepth = 100;
                } else if (newValue === 'DNA Pooling') {
                    this.sequencingDepth = 50;
                }
            }
        },
        '$data': {
            handler(newValue) {
                if (this.selectedStrategy) {
                    let nrOfSamplesForExtraction;
                    if (this.selectedStrategy === 'Soil pooling') {
                        nrOfSamplesForExtraction = this.numSites * this.numSamples / this.poolingNumber;
                    } else {
                        nrOfSamplesForExtraction = this.numSites * this.numSamples;
                    }
                    this.storageOfPlatform = nrOfSamplesForExtraction * (this.sampleSequencingDepth);
                }
            },
            deep: true
        }
    },
    mounted() {
        console.log('app mounted successfully!');
    }
}); 
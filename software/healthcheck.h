bool inmHealthCheckStart(uint8_t N, double K, bool debug);
void inmHealthCheckStop(void);
bool inmHealthCheckAddBit(bool bit, bool even);
bool inmHealthCheckOkToUseData(void);
double inmHealthCheckEstimateK(void);
double inmHealthCheckEstimateEntropyPerBit(void);

// Returns number of bytes of entropy added so far
uint32_t inmGetEntropyLevel(void);
void inmClearEntropyLevel(void);
void inmWriteEntropyToPool(uint8_t bytes, uint32_t length, uint32_t entropy);

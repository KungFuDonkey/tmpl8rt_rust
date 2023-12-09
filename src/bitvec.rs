

pub struct BitVector
{
    pub data: Vec<u64>
}

impl BitVector
{
    pub fn new(capacity: usize) -> Self
    {
        let needed_u64s = ((capacity / 64) + 1);
        BitVector
        {
            data: vec![0u64; needed_u64s]
        }
    }

    pub fn get_unchecked(&self, index: usize) -> bool
    {
        let bin = index / 64;
        let bit_index = index & 63;
        let bit_mask = 1 << bit_index;

        return self.data[bin] & bit_mask > 0;
    }

    pub fn set_true(&mut self, index: usize)
    {
        let bin = index / 64;
        let bit_index = index & 63;
        let bit_mask = 1 << bit_index;
        self.data[bin] |= bit_mask;
    }

    pub fn set_false(&mut self, index: usize)
    {
        let bin = index / 64;
        let bit_index = index & 63;
        let bit_mask = 1 << bit_index;
        let inverted_mask = u64::MAX ^ bit_mask;
        self.data[bin] &= inverted_mask;
    }

    pub fn set_unchecked(&mut self, index: usize, value: bool)
    {
        let bin = index / 64;
        let bit_index = index & 63;
        let bit_mask = 1 << bit_index;

        if value
        {
            self.data[bin] |= bit_mask;
            return;
        }

        let inverted_mask = u64::MAX ^ bit_mask;
        self.data[bin] &= inverted_mask;
    }
}
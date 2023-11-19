
pub struct Input
{
    has_focus: bool,
    key_state: [u32; 256]
}

impl Input
{
    pub fn new() -> Self
    {
        Input{
            has_focus: true,
            key_state: [0; 256]
        }
    }

    pub fn is_key_down(&self, key: u32) -> bool
    {
        return self.key_state[(key & 255) as usize] == 1;
    }

    pub fn window_has_focus(&self) -> bool
    {
        return self.has_focus;
    }

    pub fn set_key(&mut self, key: u32, value: u32)
    {
        self.key_state[(key & 255) as usize] = value;
    }

    pub fn set_focus(&mut self, focus: bool)
    {
        self.has_focus = focus;
    }

}
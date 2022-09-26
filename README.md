# Party People

Player/NPC formation movement tests

## To do

- [x] Build some sort of matrix stack system for quickly rendering all the pieces of a person
- [x] Correctly transform lighting into local space (do I need a full inverse transform? Can I maintain it in the matrix stack?)
- [X] Follow/orbit camera
- [ ] Get some level geometry going
- [ ] Keyboard control (with mouse for looking?)
- [ ] Gamepad control
- [ ] Character animation

## Notes

- Environment can be aligned with grid; can make the walls short for simplicity, or shorten the walls when the camera looks through them.
- Keep it all on a plane for simplicity

## Ideas

- Lead person A to person B
- Lead person A to destination without person A being seen
- Guard a door
- Forced-march to a destination
- Gang isolating and surrounding player
- Buddy movement: single-file, side-by-side, back-to-back
- Transition movement types based on whether there are doors or enemies
- Everyone follows player
- NPC has constraints on where they want to go

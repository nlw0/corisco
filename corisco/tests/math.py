from attest import TestBase, test, Assert

class Math(TestBase):

    def __context__(self):
        self.value = 1 + 1
        yield

    @test
    def arithmetics(self):
        Assert(self.value) == 3

if __name__ == '__main__':
    from attest import Tests
    suite = Tests([Math()])
    suite.run()

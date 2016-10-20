from __future__ import print_function
import unittest
import datetime
from transcriptic_tools.custom_connection import CustomConnection as Connection

def dict_subset(d, keys):
    return dict((k, d[k]) for k in keys)

class TestCustomConnection(unittest.TestCase):
    
    def setUp(self):
        self.api = Connection.from_file(".transcriptic")
    
    def test_run(self):
        run = self.api.run('p18v89yurexft', 'r195cxaehbw3m')
        
        self.assertEqual(dict_subset(run,['id','status','title','created_at','completed_at',
                                          'conversation_id']),
                         {
                             'status': 'complete',
                             'title': 'UnitTest Run (Do Not Edit)',
                             'created_at': '2016-06-22T14:31:44.925-07:00',
                             'completed_at': '2016-06-22T14:33:00.670-07:00',
                             'conversation_id': 'conv195cxaep8pk3',
                             'id': 'r195cxaehbw3m'
                         }
                         )      
        
    def test_create_post(self):
        """
        To manually delete posts if there are issues, go here
        https://secure.transcriptic.com/becker-lab/p18v89yurexft/runs/r195cxaehbw3m
        
        """
        
        curr_time = datetime.datetime.now().strftime("%m_%d_%Y_%H_%M_%S")
        
        post_text = "Test Comment %s"%curr_time
        
        project_id = 'p18v89yurexft'
        run_id = 'r195cxaehbw3m'
        
        post = self.api.create_post(project_id, run_id, post_text)
        
        posts = self.api.posts(project_id, run_id)
        
        self.assertEqual(len(posts),1)
        self.assertEqual(posts[0]['text'], post_text)
        
        self.api.delete_post(project_id,run_id,post['id'])
        
        posts_after_delete = self.api.posts(project_id, run_id)
        
        self.assertEqual(posts_after_delete,[])
        
        
    
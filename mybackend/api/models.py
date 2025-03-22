from django.db import models
from django.contrib.auth.models import AbstractUser


class CustomUser(AbstractUser):
    pass  # 可以添加额外字段，如手机号等

class QueryHistory(models.Model):
    user = models.ForeignKey(CustomUser, on_delete=models.CASCADE)  # 关联用户
    smiles = models.CharField(max_length=255)  # 存储 SMILES 字符串
    created_at = models.DateTimeField(auto_now_add=True)  # 记录查询时间

    def __str__(self):
        return f"{self.user.username} - {self.smiles}"


